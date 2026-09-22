#' Initialize a group of server modules on its first visible tab.
#'
#' The callback runs in the session domain and is isolated from later input
#' changes. Modules retain their observers and state when the user leaves a tab.
#' @noRd
on_first_tab <- function(input, tabs, initialize) {
  force(tabs)
  force(initialize)
  observer <- shiny::observeEvent(input$navbarID, {
    if (!input$navbarID %in% tabs) return()
    initialize()
    observer$destroy()
  }, ignoreInit = FALSE, priority = -1)
  invisible(observer)
}
