check_plot_on_start <- function(browser_options) {
  !as.logical(browser_options["plot_on_start"]) ||
    (!is.null(isolate(getQueryString())[["go"]]) && isolate(getQueryString())[["go"]] == TRUE)
}
