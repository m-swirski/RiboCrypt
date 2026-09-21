check_plot_on_start <- function(browser_options) {
  url_go <- isolate(getQueryString())[["go"]]
  !isTRUE(as.logical(browser_options["plot_on_start"])) ||
    isTRUE(as.logical(url_go))
}
