#' Reuse the annotation track's exact transcript-coordinate intervals.
#' @noRd
megabrowser_translon_panel <- function(controller) {
  annotation <- annotation_controller(controller$dff, controller$display_region,
    controller$annotation, leader_extension = 0, trailer_extension = 0, viewMode = controller$viewMode)
  createGeneModelPanel(annotation$display_range, annotation$annotation,
    tx_annotation = controller$tx_annotation, custom_regions = controller$customRegions,
    viewMode = controller$viewMode, collapse_intron_flank = controller$collapsed_introns_width)[[1]]
}

#' Deduplicate complete exon sets, not just their outer bounds.
#' @noRd
megabrowser_translon_regions <- function(panel, cds_names) {
  coding <- panel[type == "cds"]
  validate(need(nrow(coding) > 0L, "No coding regions are visible in the annotation track."))
  sets <- lapply(split(coding, coding$gene_names), function(rows) IRanges::reduce(IRanges::IRanges(rows$rect_starts, rows$rect_ends)))
  keys <- vapply(sets, function(ranges) paste(start(ranges), end(ranges), sep = ":", collapse = ";"), character(1))
  regions <- lapply(unique(keys), function(key) {
    labels <- names(sets)[keys == key]
    list(labels = paste(labels, collapse = ", "), kind = if (any(labels %in% cds_names)) "CDS" else "Translon",
         ranges = sets[[match(key, keys)]])
  })
  megabrowser_clean_cds(regions)
}

#' Remove upstream ORF overlap from each CDS without filling exon gaps.
#' @noRd
megabrowser_clean_cds <- function(regions) {
  cds <- which(vapply(regions, function(region) region$kind == "CDS", logical(1)))
  for (index in cds) {
    original <- regions[[index]]
    upstream <- Filter(function(region) region$kind == "Translon" && min(start(region$ranges)) < min(start(original$ranges)), regions)
    if (!length(upstream)) next
    overlap <- IRanges::reduce(do.call(c, lapply(upstream, `[[`, "ranges")))
    clean <- IRanges::setdiff(original$ranges, overlap)
    if (!identical(clean, original$ranges)) regions[[length(regions) + 1L]] <-
      list(labels = paste0("clean_cds (", original$labels, ")"), kind = "clean_cds", ranges = clean,
           parent_labels = original$labels, parent_start = min(start(original$ranges)), parent_end = max(end(original$ranges)))
  }
  names(regions) <- paste0("R", seq_along(regions))
  regions
}

#' Region means on unbinned raw coverage; normalized heatmap values are not ratios.
#' @noRd
megabrowser_translon_density <- function(raw, regions) {
  positions <- seq_len(nrow(raw))
  values <- vapply(regions, function(region) {
    keep <- IRanges::overlapsAny(IRanges::IRanges(positions, width = 1), region$ranges)
    if (!any(keep)) return(rep(NA_real_, ncol(raw)))
    colMeans(raw[keep, , drop = FALSE])
  }, numeric(ncol(raw)))
  matrix(values, nrow = ncol(raw), dimnames = list(colnames(raw), names(regions)))
}

#' Canonical pair orientation keeps CDS references in the denominator.
#' @noRd
megabrowser_translon_pairs <- function(regions) {
  validate(need(length(regions) >= 2L, "At least two distinct coding regions are needed. Enable predicted translons and Generate Plot."))
  pairs <- t(utils::combn(names(regions), 2L))
  priority <- vapply(regions, function(region) match(region$kind, c("Translon", "User defined", "CDS", "clean_cds")), integer(1))
  swap <- priority[pairs[, 1]] > priority[pairs[, 2]]
  pairs[swap, ] <- pairs[swap, 2:1, drop = FALSE]
  labels <- vapply(regions, `[[`, character(1), "labels")
  data.table(Pair = paste(pairs[, 1], pairs[, 2], sep = "/"), Numerator = pairs[, 1], Denominator = pairs[, 2],
    Numerator_label = unname(labels[pairs[, 1]]), Denominator_label = unname(labels[pairs[, 2]]))
}

#' Finite ratios only: zero denominators remain undefined, without pseudocounts.
#' @noRd
megabrowser_translon_ratios <- function(density, pairs, groups) {
  membership <- rep(names(groups), lengths(groups))
  names(membership) <- rownames(density)[unlist(groups, use.names = FALSE)]
  rbindlist(lapply(seq_len(nrow(pairs)), function(index) {
    numerator <- density[, pairs$Numerator[index]]
    denominator <- density[, pairs$Denominator[index]]
    ratio <- ifelse(is.finite(denominator) & denominator > 0 & is.finite(numerator), numerator / denominator, NA_real_)
    data.table(Pair = pairs$Pair[index], Run = rownames(density), Group = unname(membership[rownames(density)]),
      Numerator_label = pairs$Numerator_label[index], Denominator_label = pairs$Denominator_label[index],
      Numerator_density = numerator, Denominator_density = denominator, Ratio = ratio,
      Log2_ratio = ifelse(is.finite(ratio) & ratio > 0, log2(ratio), NA_real_))
  }))
}

#' Exploratory cluster-versus-rest effect sizes; libraries are not biological replicates.
#' @noRd
megabrowser_translon_test <- function(values, rest) {
  if (length(values) < 2L || length(rest) < 2L || length(unique(c(values, rest))) < 2L)
    return(list(P = NA_real_, Rank_biserial = NA_real_))
  test <- suppressWarnings(stats::wilcox.test(values, rest, exact = FALSE))
  list(P = test$p.value, Rank_biserial = 2 * unname(test$statistic) / (length(values) * length(rest)) - 1)
}

#' Summaries retain zeros and explicitly count undefined ratios.
#' @noRd
megabrowser_translon_stats <- function(ratios) {
  result <- ratios[, {
    values <- Ratio[is.finite(Ratio)]
    rest <- ratios[Pair == .BY$Pair & Group != .BY$Group & is.finite(Ratio), Ratio]
    quantiles <- if (length(values)) unname(quantile(values, c(0.25, 0.5, 0.75))) else rep(NA_real_, 3)
    c(list(Libraries = .N, Valid = length(values), Undefined = sum(!is.finite(Ratio)),
           Zero = sum(values == 0), Q25 = quantiles[1], Median = quantiles[2], Q75 = quantiles[3],
           `Ratio > 1 (%)` = if (length(values)) 100 * mean(values > 1) else NA_real_),
      megabrowser_translon_test(values, rest))
  }, by = .(Pair, Group)]
  result[, BH := p.adjust(P, method = "BH")]
  result
}

#' Small region table includes aliases, coordinates and removed CDS bases.
#' @noRd
megabrowser_translon_region_table <- function(regions) {
  rbindlist(lapply(names(regions), function(id) {
    region <- regions[[id]]
    data.table(Region = id, Type = region$kind, Labels = region$labels,
      Length_nt = sum(width(region$ranges)), Intervals = paste(start(region$ranges), end(region$ranges), sep = "-", collapse = ";"))
  }))
}

#' Analyze original library membership and raw coverage only on demand.
#' @noRd
megabrowser_translon_analysis <- function(controller, table, grouped, custom = list(), raw = NULL) {
  regions <- megabrowser_translon_regions(megabrowser_translon_panel(controller), names(controller$annotation))
  regions <- c(regions, custom)
  names(regions) <- paste0("R", seq_along(regions))
  pairs <- megabrowser_translon_pairs(regions)
  if (is.null(raw)) raw <- as.matrix(load_collection(controller$table_path, columns = colnames(table$table)))
  raw <- raw[, colnames(table$table), drop = FALSE]
  validate(need(all(is.finite(raw)) && all(raw >= 0), "Translon ratios require finite, non-negative raw coverage."))
  bounds <- unlist(lapply(regions, function(region) c(start(region$ranges), end(region$ranges))), use.names = FALSE)
  validate(need(all(bounds >= 1 & bounds <= nrow(raw)), "Annotation coordinates do not match the raw coverage range."))
  density <- megabrowser_translon_density(raw, regions)
  ratios <- megabrowser_translon_ratios(density, pairs, megabrowser_display_groups(grouped$meta, colnames(raw)))
  statistics <- merge(megabrowser_translon_stats(ratios), pairs[, .(Pair, Numerator_label, Denominator_label)], by = "Pair", sort = FALSE)
  list(regions = megabrowser_translon_region_table(regions), pairs = pairs, ratios = ratios, statistics = statistics)
}

#' Strict one-based inclusive coordinates for an additional transcript region.
#' @noRd
megabrowser_user_region <- function(coordinates, label, limit) {
  coordinates <- trimws(coordinates)
  if (length(coordinates) != 1L || is.na(coordinates) || !grepl("^[0-9]+(\\s*:\\s*[0-9]+)?$", coordinates))
    stop("Enter a position or start:end, for example 40 or 40:80.", call. = FALSE)
  bounds <- as.numeric(strsplit(coordinates, ":", fixed = TRUE)[[1]])
  if (length(bounds) == 1L) bounds <- rep(bounds, 2L)
  if (any(!is.finite(bounds)) || any(bounds < 1 | bounds > limit))
    stop(paste0("Coordinates must be between 1 and ", limit, "."), call. = FALSE)
  if (bounds[1] > bounds[2]) stop("Start must not exceed end.", call. = FALSE)
  label <- trimws(label)
  if (length(label) != 1L || is.na(label) || !nzchar(label)) stop("Enter a region label.", call. = FALSE)
  list(labels = label, kind = "User defined", ranges = IRanges::IRanges(bounds[1], bounds[2]))
}
