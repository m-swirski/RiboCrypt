#' Export the exact profiles retained with the generated plot, not current inputs.
#' @noRd
browser_coverage_table <- function(plot) {
  snapshot <- attr(plot, "coverage_export")
  shiny::req(length(snapshot$profiles) > 0)
  profiles <- snapshot$profiles
  positions <- profiles[[1]]$position
  stopifnot(all(vapply(profiles, function(x) identical(x$position, positions), logical(1))))
  values <- lapply(profiles, function(x) x$count)
  labels <- snapshot$labels
  if (isTRUE(snapshot$summary)) {
    values <- c(list(Reduce(`+`, values)), rev(values))
    labels <- c("summary", rev(labels))
  }
  data.table::as.data.table(stats::setNames(c(list(positions), values), c("position", labels)))
}

#' @noRd
browser_coverage_filename <- function(controls) {
  tx <- names(controls$display_region)[1]
  if (!length(tx) || is.na(tx) || !nzchar(tx)) tx <- "region"
  paste0("RiboCrypt_", gsub("[^A-Za-z0-9_.-]", "_", tx), "_coverage.csv")
}

#' Register lazily: no coverage extraction or serialization during startup.
#' @noRd
browser_coverage_download <- function(output, controls, browser_plot) {
  output$download_coverage <- rc_download_handler(
    filename = function() browser_coverage_filename(controls()),
    contentType = "text/csv",
    content = function(file) {
      data.table::fwrite(browser_coverage_table(browser_plot()), file)
    }
  )
}
