#' Stable bins, including constant or missing numeric metadata.
#' @noRd
megabrowser_numeric_bins <- function(values, bins = 5L) {
  finite <- is.finite(values)
  result <- rep("(Missing)", length(values))
  if (any(finite)) result[finite] <- if (length(unique(values[finite])) == 1L) "All" else
    as.character(cut(values[finite], breaks = bins, include.lowest = TRUE))
  factor(result, levels = unique(result[order(values, na.last = TRUE)]))
}

#' Metadata in the same library order as the clustered full matrix.
#' @noRd
megabrowser_sample_metadata <- function(table, enrichment_test_on, numeric_bins = 5L) {
  values <- table$metadata_field
  orders <- unlist(attr(table$table, "row_order_list"), use.names = FALSE)
  meta <- data.table(Run = colnames(table$table), grouping = unname(values), order = seq_along(values))
  other <- attr(values, "other_columns")
  for (column in setdiff(names(other), names(meta))) meta[[column]] <- other[[column]][match(meta$Run, other$Run)]
  meta <- meta[orders]
  row_orders <- attr(table$table, "row_order_list")
  if (length(row_orders) > 1L) meta[, cluster := rep(seq_along(row_orders), lengths(row_orders))]
  megabrowser_metadata_groups(meta, values, enrichment_test_on, numeric_bins)
}

#' Attach analysis groups without reducing the number of libraries.
#' @noRd
megabrowser_metadata_groups <- function(meta, values, term, bins) {
  numeric_group <- is.numeric(meta$grouping)
  if (term %in% c("Ratio bins", "Other gene tpm bins")) {
    meta[, cluster := megabrowser_numeric_bins(grouping, bins)]
    other <- attr(values, "other_columns")
    field <- setdiff(names(other), "Run")[1L]
    meta[, grouping := other[[field]][match(Run, other$Run)]]
  } else if (numeric_group) meta[, grouping_numeric_bins := megabrowser_numeric_bins(grouping, bins)]
  label <- attr(values, "xlab")
  suffix <- if (numeric_group && !term %in% c("Ratio bins", "Other gene tpm bins")) " (Numeric Bins)" else
    if (!is.null(meta$cluster) && !is.factor(meta$cluster)) " Clusters (K-means)" else ""
  attr(meta, "xlab") <- paste0(label, suffix)
  attr(meta, "ylab") <- if (is.null(meta$cluster)) "Counts" else "Enrichment"
  attr(meta, "runIDs") <- attr(values, "runIDs")[match(meta$Run, attr(values, "runIDs")$Run)]
  meta
}

#' Recompute enrichment for one metadata field using the full sample membership.
#' @noRd
megabrowser_metadata_enrichment <- function(grouped, metadata, field = "grouping") {
  if (is.null(field) || identical(field, "") || identical(field, "grouping")) {
    if (is.null(grouped$enrich_dt)) grouped$enrich_dt <- allsamples_meta_stats(grouped$meta)
    return(grouped)
  }
  stopifnot(field %in% names(metadata))
  meta <- copy(grouped$meta)
  values <- metadata[[field]][match(meta$Run, metadata$Run)]
  if ("grouping_numeric_bins" %in% names(meta)) meta[, grouping_numeric_bins := NULL]
  meta[, grouping := if (is.numeric(values)) megabrowser_numeric_bins(values) else
    fifelse(is.na(values) | as.character(values) == "", "(Missing)", as.character(values))]
  attr(meta, "xlab") <- field
  list(meta = meta, enrich_dt = allsamples_meta_stats(meta))
}

#' Display groups follow the existing clusters, numeric bins or ordering categories.
#' @noRd
megabrowser_display_groups <- function(meta, runs) {
  column <- if ("cluster" %in% names(meta)) "cluster" else
    if ("grouping_numeric_bins" %in% names(meta)) "grouping_numeric_bins" else "grouping"
  groups <- as.character(meta[[column]])
  groups[is.na(groups) | groups == ""] <- "(Missing)"
  ordered <- groups[match(runs, meta$Run)]
  split(seq_along(runs), factor(ordered, levels = rev(unique(groups))))
}

#' Retain complete groups for display without changing sample-level analysis.
#' @noRd
megabrowser_focused_display <- function(table, meta, selected = NULL) {
  groups <- megabrowser_display_groups(meta, colnames(table))
  selected <- intersect(names(groups), selected)
  if (!length(selected) || length(selected) == length(groups)) return(list(table = table, meta = meta))
  keep <- sort(unlist(groups[selected], use.names = FALSE))
  focused <- subset_collection_columns(table, colnames(table)[keep])
  focused <- megabrowser_focus_attributes(focused, table, keep)
  sidebar <- copy(meta[Run %in% colnames(focused)])
  sidebar[, index := .I]
  list(table = focused, meta = sidebar)
}

#' Remap existing display orders rather than reclustering a focused matrix.
#' @noRd
megabrowser_focus_attributes <- function(focused, original, keep) {
  orders <- lapply(attr(original, "row_order_list"), function(order) match(intersect(order, keep), keep))
  orders <- orders[lengths(orders) > 0L]
  km <- attr(original, "km")
  km$cluster <- km$cluster[keep]
  attr(focused, "km") <- km
  attr(focused, "row_order_list") <- orders
  attr(focused, "clusters") <- length(orders)
  focused
}

#' A compact metadata value for the sidebar of a collapsed group.
#' @noRd
megabrowser_group_value <- function(values) {
  if (is.numeric(values)) return(if (all(is.na(values))) NA_real_ else mean(values, na.rm = TRUE))
  values <- as.character(values)
  values <- values[!is.na(values) & nzchar(values)]
  counts <- sort(table(values, useNA = "no"), decreasing = TRUE)
  if (!length(counts)) return("(Missing)")
  names(counts)[1L]
}

#' Collapse only the display matrix; original counts, clustering and statistics stay intact.
#' @noRd
megabrowser_collapsed_display <- function(table, meta) {
  groups <- megabrowser_display_groups(meta, colnames(table))
  matrix <- if (is.matrix(table)) table else as.matrix(table)
  collapsed <- vapply(groups, function(indices) matrixStats::rowMeans2(matrix, cols = indices), numeric(nrow(matrix)))
  collapsed <- matrix(collapsed, nrow = nrow(matrix), dimnames = list(rownames(matrix), names(groups)))
  collapsed <- set_collection_user_attributes(collapsed, collection_user_attributes(table))
  attr(collapsed, "km") <- list(cluster = seq_along(groups))
  attr(collapsed, "row_order_list") <- as.list(seq_along(groups))
  attr(collapsed, "clusters") <- length(groups)
  attr(collapsed, "collapsed_clusters") <- TRUE
  list(table = collapsed, meta = megabrowser_collapsed_sidebar(meta, colnames(matrix), groups))
}

#' Sidebar rows use the same group order and member counts as the collapsed matrix.
#' @noRd
megabrowser_collapsed_sidebar <- function(meta, runs, groups) {
  columns <- setdiff(names(meta), c("Run", "order", "index", "cluster"))
  rows <- lapply(groups, function(indices) {
    members <- meta[match(runs[indices], Run)]
    as.data.table(lapply(members[, ..columns], megabrowser_group_value))
  })
  result <- rbindlist(rows)
  result[, `:=`(cluster = names(groups), libraries = lengths(groups))]
  result <- result[rev(seq_len(.N))]
  result[, index := .I]
  attr(result, "xlab") <- attr(meta, "xlab")
  result
}
