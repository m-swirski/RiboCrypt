#' Original library indices behind each displayed row.
#' @noRd
megabrowser_cell_memberships <- function(view, display, meta) {
  if (isTRUE(attr(display$table, "collapsed_clusters")))
    return(megabrowser_display_groups(meta, colnames(view))[colnames(display$table)])
  setNames(lapply(colnames(display$table), function(run) match(run, colnames(view))), colnames(display$table))
}

#' Descriptive support gate uses coverage sums, not independent read counts.
#' @noRd
megabrowser_region_support <- function(view, display, meta, minimum = 10, libraries = 3) {
  validate(need(length(minimum) == 1L && is.finite(minimum) && minimum >= 0, "Minimum coverage sum must be a non-negative number."))
  validate(need(length(libraries) == 1L && is.finite(libraries) && libraries >= 1 && libraries == floor(libraries), "Minimum supporting libraries must be a positive integer."))
  lengths <- attr(view, "translon_regions")$Length_nt
  sums <- sweep(view, 1L, lengths, "*")
  adequate <- is.finite(sums) & sums > 0 & sums >= minimum
  groups <- megabrowser_cell_memberships(view, display, meta)
  support <- vapply(groups, function(indices) rowSums(adequate[, indices, drop = FALSE]), numeric(nrow(view)))
  support <- matrix(support, nrow = nrow(view), dimnames = list(rownames(view), names(groups)))
  list(count = support, low = support < if (isTRUE(attr(display$table, "collapsed_clusters"))) libraries else 1L)
}

#' Overlay only glyphs; heatmap values, colours and saved clustering are untouched.
#' @noRd
megabrowser_support_plot <- function(plot, view, display, meta, input) {
  if (!isTRUE(attr(display$table, "collapsed_translons")) || !isTRUE(input$region_support)) return(plot)
  support <- megabrowser_region_support(view, display, meta, input$region_min_coverage, input$region_min_libraries)
  orders <- unlist(attr(display$table, "row_order_list"), use.names = FALSE)
  points <- which(t(support$low)[orders, , drop = FALSE], arr.ind = TRUE)
  if (!nrow(points)) return(plot)
  plot <- plotly::plotly_build(plot)
  plot$x$data[[length(plot$x$data) + 1L]] <- list(x = points[, "col"], y = points[, "row"],
    type = if (nrow(points) > 2000) "scattergl" else "scatter", mode = "markers",
    marker = list(symbol = "x", size = 6, color = "#666666"), hoverinfo = "skip",
    showlegend = FALSE, xaxis = "x", yaxis = "y")
  plot
}
