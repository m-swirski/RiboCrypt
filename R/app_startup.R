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

#' Load the normal Plotly dependencies before an automatic startup plot arrives.
#' @noRd
plotly_startup_dependencies <- function(plot_on_start) {
  if (!isTRUE(as.logical(plot_on_start))) return(NULL)
  plotly::plot_ly()$dependencies
}

# Prime Plotly for Observatory. Eager loading favors the first Observatory
# transition; prefetch avoids blocking the initial page, at a speed cost.
rc_observatory_plotly_preload <- function(enabled, plot_on_start) {
  if (!isTRUE(enabled) || isTRUE(as.logical(plot_on_start))) return(NULL)
  mode <- getOption("ribocrypt.observatory_plotly_preload", "eager")
  if (identical(mode, "none")) return(NULL)
  if (identical(mode, "eager")) return(plotly::plot_ly()$dependencies)
  if (!identical(mode, "prefetch")) stop("Invalid ribocrypt.observatory_plotly_preload option")
  deps <- plotly::plot_ly()$dependencies
  match <- Filter(function(x) identical(x$name, "plotly-main"), deps)
  if (length(match) != 1L) return(NULL)
  dep <- match[[1L]]
  if (length(dep$script) != 1L || is.null(dep$package)) return(NULL)
  dir <- system.file(dep$src$file, package=dep$package)
  if (!nzchar(dir) || !file.exists(file.path(dir, dep$script))) return(NULL)
  prefix <- paste0(dep$name, "-", as.character(dep$version))
  shiny::addResourcePath(prefix, dir)
  shiny::tags$link(rel="prefetch", as="script", href=paste0(prefix,"/",dep$script))
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
