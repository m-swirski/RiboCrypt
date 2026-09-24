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

#' Cache the rendered shell of this app's request-independent, static UI.
#'
#' Keep Shiny's complete response, including dependency/singleton registration.
#' Query requests (including bookmarks/showcase) and non-UI requests bypass it.
#' This must not be used for a request-dependent UI function.
#' @noRd
cache_static_app_ui <- function(app) {
  handler <- app$httpHandler
  cached <- NULL
  app$httpHandler <- function(req) {
    plain <- identical(req$REQUEST_METHOD, "GET") && identical(req$PATH_INFO, "/") &&
      !nzchar(req$QUERY_STRING %||% "")
    if (!plain) return(handler(req))
    if (!is.null(cached)) return(cached)
    response <- handler(req)
    if (!is.null(response) && isTRUE(response$status == 200L)) cached <<- response
    response
  }
  app
}
