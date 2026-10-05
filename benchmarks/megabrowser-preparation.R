# Source after devtools::load_all("."). Synthetic backend benchmark; not app latency.
local({
  baseline_data <- function(table) {
    mat <- RiboCrypt:::megabrowser_ordered_matrix(table)
    ratio <- RiboCrypt:::megabrowser_table_ratio(table)
    list(x = ((seq_len(ncol(mat)) - 1L) * ratio) + 1L, y = seq_len(nrow(mat)),
      type = RiboCrypt:::megabrowser_heatmap_renderer(length(mat), isTRUE(attr(table, "collapsed_clusters"))),
      z = mat[rev(seq_len(nrow(mat))), , drop = FALSE],
      x_range = RiboCrypt:::megabrowser_full_x_range(table = table), y_range = c(0.5, nrow(mat) + 0.5))
  }
  measure <- function(fun, table, repetitions = 10L) {
    fun(table)
    vapply(seq_len(5), function(trial) {
      gc()
      unname(system.time(for (i in seq_len(repetitions)) fun(table))["elapsed"]) / repetitions
    }, numeric(1))
  }
  results <- lapply(list(c(483L, 1754L), c(2000L, 3000L)), function(size) {
    set.seed(42)
    table <- matrix(runif(prod(size)), size[1], dimnames = list(NULL, paste0("r", seq_len(size[2]))))
    attr(table, "row_order_list") <- split(sample(seq_len(size[2])), rep(1:5, length.out = size[2]))
    attr(table, "km") <- list(cluster = rep(1:5, length.out = size[2]))
    stopifnot(identical(baseline_data(table), RiboCrypt:::megabrowser_plotly_heatmap_data(table)))
    list(positions = size[1], libraries = size[2],
      before = measure(baseline_data, table), after = measure(RiboCrypt:::megabrowser_plotly_heatmap_data, table))
  })
  jsonlite::write_json(list(R = R.version.string, results = results),
    "/tmp/megabrowser-preparation.json", auto_unbox = TRUE, pretty = TRUE)
  print(results)
})
