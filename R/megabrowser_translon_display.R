#' Equal-width region columns, computed before any library-group collapse.
#' @noRd
megabrowser_translon_display <- function(table, raw, regions) {
  regions <- megabrowser_ordered_translon_regions(regions)
  validate(need(length(regions) > 0L, "No translon regions are available for display."))
  raw <- raw[, colnames(table), drop = FALSE]
  bounds <- unlist(lapply(regions, function(region) c(start(region$ranges), end(region$ranges))), use.names = FALSE)
  validate(need(all(bounds >= 1 & bounds <= nrow(raw)), "Region coordinates do not match raw coverage."))
  validate(need(all(is.finite(raw)) && all(raw >= 0), "Region coverage must be finite and non-negative."))
  split_frames <- isTRUE(attr(table, "split_by_frame"))
  if (split_frames) {
    # Work from unbinned raw coverage: each phase uses only its own bases,
    # including discontiguous regions and regions shorter than a codon.
    positions <- seq_len(nrow(raw))
    position_ranges <- IRanges::IRanges(positions, width = 1)
    region_info <- megabrowser_translon_region_table(regions)
    region_info$Display_label <- megabrowser_biological_labels(regions)
    values <- vector("list", length(regions) * 3L)
    info <- vector("list", length(values))
    index <- 0L
    for (i in seq_along(regions)) {
      selected <- IRanges::overlapsAny(position_ranges, regions[[i]]$ranges)
      # A clean CDS can start after an excluded overlap; retain its parent's
      # coding frame rather than redefining frame zero at the first kept base.
      origin <- if (!is.null(regions[[i]]$parent_start)) regions[[i]]$parent_start else
        if (length(regions[[i]]$ranges)) min(start(regions[[i]]$ranges)) else NA_integer_
      phases <- (positions - origin) %% 3L
      for (frame in 0:2) {
        index <- index + 1L
        keep <- selected & !is.na(phases) & phases == frame
        values[[index]] <- if (any(keep)) colMeans(raw[keep, , drop = FALSE]) else rep(NA_real_, ncol(raw))
        row <- data.table::copy(region_info[i, ])
        row[, `:=`(Parent_region = Region, Region = paste0(Region, "_F", frame),
                    Frame = frame, Frame_origin = origin, Length_nt = sum(keep),
                    Labels = paste0(Labels, " (frame ", frame, ")"))]
        info[[index]] <- row
      }
    }
    density <- do.call(rbind, values)
    region_info <- data.table::rbindlist(info)
    rownames(density) <- region_info$Region
    colnames(density) <- colnames(raw)
  } else {
    density <- t(megabrowser_translon_density(raw, regions))
    region_info <- megabrowser_translon_region_table(regions)
    region_info$Display_label <- megabrowser_biological_labels(regions)
  }
  density <- set_collection_user_attributes(density, collection_user_attributes(table))
  attr(density, "ratio") <- 1L
  attr(density, "collapsed_translons") <- TRUE
  attr(density, "translon_regions") <- region_info
  attr(density, "translon_summary") <- rowMeans(density)
  density
}

#' Region density summary uses individual libraries, not equally weighted groups.
#' @noRd
megabrowser_translon_summary_plot <- function(table) {
  regions <- attr(table, "translon_regions")
  bar_color <- if ("Frame" %in% names(regions)) c("#e3524a", "#3aa35b", "#3f7fd6")[regions$Frame + 1L] else "#198c83"
  plot <- plotly::plot_ly() %>% plotly::add_bars(x = seq_len(nrow(table)), y = attr(table, "translon_summary"),
    width = 1, marker = list(color = bar_color), textposition = "none",
    text = regions$Labels,
    hovertemplate = "%{text}<br>Mean coverage: %{y}<extra></extra>")
  mb_finalize_summary_track(plot, megabrowser_full_x_range(table = table)) %>%
    plotly::layout(margin = megabrowser_translon_margins())
}

#' Region labels replace genomic geometry when each region has equal width.
#' @noRd
megabrowser_translon_annotation_plot <- function(table) {
  regions <- attr(table, "translon_regions")
  labels <- regions$Display_label
  if ("Frame" %in% names(regions)) labels[regions$Frame != 1L] <- ""
  plot <- plotly::plot_ly() %>% plotly::add_bars(x = seq_len(nrow(table)), y = rep(1, nrow(table)), width = 1,
    marker = list(color = ifelse(regions$Type %in% c("CDS", "clean_cds"), "#2878a5", "#198c83")),
    text = labels, textposition = "inside", customdata = paste(regions$Region, regions$Labels, regions$Intervals),
    hovertemplate = "%{text}: %{customdata}<extra></extra>") %>%
    plotly::layout(margin = megabrowser_translon_margins(), bargap = 0,
      xaxis = mb_summary_xaxis(megabrowser_full_x_range(table = table)),
      yaxis = list(visible = FALSE, range = c(0, 1), fixedrange = TRUE), showlegend = FALSE) %>%
    plotly::config(doubleClick = FALSE)
  plot$x$source <- "mb_bottom"
  plot
}

#' Coordinate-based names do not imply independent translation evidence.
#' @noRd
megabrowser_biological_labels <- function(regions) {
  cds <- Filter(function(r) r$kind %in% c("CDS", "clean_cds"), regions)
  first <- if (length(cds)) min(vapply(cds, function(r) if (!is.null(r$parent_start)) r$parent_start else min(start(r$ranges)), numeric(1))) else NA_real_
  last <- if (length(cds)) max(vapply(cds, function(r) if (!is.null(r$parent_end)) r$parent_end else max(end(r$ranges)), numeric(1))) else NA_real_
  labels <- vapply(regions, function(r) {
    if (r$kind == "User defined") return(r$labels)
    if (r$kind == "clean_cds") return("clean CDS")
    if (r$kind == "CDS") return("CDS")
    if (is.na(first)) return("ORF")
    if (min(start(r$ranges)) < first) return("uORF")
    if (min(start(r$ranges)) > last) "dORF" else "internal ORF"
  }, character(1))
  for (kind in intersect(c("uORF", "dORF", "internal ORF", "ORF"), labels)) labels[labels == kind] <- paste0(kind, seq_len(sum(labels == kind)))
  unname(labels)
}

#' Reserve a shared margin for library indices, without axis auto-expansion.
#' @noRd
megabrowser_translon_margins <- function() {
  margins <- megabrowser_heatmap_margins()
  margins$l <- 50
  margins$autoexpand <- FALSE
  margins
}

#' Prefer each clean CDS over its parent and retain transcript-coordinate order.
#' @noRd
megabrowser_ordered_translon_regions <- function(regions) {
  clean <- Filter(function(region) region$kind == "clean_cds", regions)
  parents <- vapply(clean, function(region) region$parent_labels, character(1))
  regions <- Filter(function(region) !(region$kind == "CDS" && region$labels %in% parents), regions)
  starts <- vapply(regions, function(region) if (length(region$ranges)) min(start(region$ranges)) else
    if (!is.null(region$parent_start)) region$parent_start else Inf, numeric(1))
  ends <- vapply(regions, function(region) if (length(region$ranges)) max(end(region$ranges)) else Inf, numeric(1))
  regions[order(starts, ends, names(regions))]
}

#' Region-specific fold change uses the full individual-library reference.
#' @noRd
megabrowser_translon_cluster_score <- function(display) {
  table <- display$table
  reference <- attr(table, "translon_summary")
  ratio <- sweep(table, 1L, reference, "/")
  ratio[!is.finite(reference) | reference <= 0, ] <- NA_real_
  score <- log2(ratio)
  attr(table, "translon_density") <- as.matrix(table)
  attr(table, "translon_fold_change") <- ratio
  table[] <- pmax(-1, pmin(1, score))
  attr(table, "translon_score") <- TRUE
  display$table <- table
  display
}

#' Compact legend uses the same palette interpolation as the heatmap.
#' @noRd
megabrowser_score_legend <- function(colors) {
  values <- c(-1, 0, 1)
  swatches <- circlize::colorRamp2(seq(-1, 1, length.out = length(colors)), colors, space = "RGB")(values)
  tagList("log2 FC", lapply(seq_along(values), function(i)
    tags$span(tags$i(style = paste0("background:", swatches[i])), c("<= 0.5x", "1x", ">= 2x")[i])))
}
