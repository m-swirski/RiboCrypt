#' Initialize a group of server modules on its first visible tab.
#'
#' The callback runs in the session domain and is isolated from later input
#' changes. Modules retain their observers and state when the user leaves a tab.
#' @noRd
on_first_tab <- function(input, tabs, initialize, prewarm = NULL,
                         session = shiny::getDefaultReactiveDomain()) {
  force(tabs)
  force(initialize)
  done <- FALSE
  run <- function() {
    if (done) return(invisible())
    done <<- TRUE
    initialize()
  }
  observer <- shiny::observeEvent(input$navbarID, {
    if (!input$navbarID %in% tabs) return()
    run()
    observer$destroy()
  }, ignoreInit = FALSE, priority = -1)
  # Optional idle prewarm: build the module servers `prewarm` seconds after the first page flush, so the first visit
  # to the tab does not pay for creating them (~0.2-0.3 s for Observatory). Hidden outputs stay suspended until shown.
  if (!is.null(prewarm) && !is.null(session)) {
    session$onFlushed(function() {
      later::later(function() {
        if (done || session$isClosed()) return()
        shiny::withReactiveDomain(session, shiny::isolate(run()))
        observer$destroy()
      }, delay = prewarm)
    }, once = TRUE)
  }
  invisible(observer)
}

#' Delay (s) after the first page flush before collection pages are prepared in the background, or NULL to disable.
#' @noRd
rc_observatory_prewarm_delay <- function() {
  delay <- getOption("ribocrypt.observatory_prewarm", 1)
  if (is.null(delay) || isFALSE(delay)) return(NULL)
  as.numeric(delay)
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
# App URL parameters (?go=FALSE, ?dff=..&gene=.., obs_state) are read in the browser and leave the static page HTML
# unchanged; Shiny's own bookmark-restore and showcase parameters change the page and always render it.
rc_ui_query_cacheable <- function(query) {
  if (!nzchar(query)) return(TRUE)
  if (!identical(shiny::getShinyOption("bookmarkStore", default = "disable"), "disable")) return(FALSE)
  keys <- sub("=.*$", "", strsplit(sub("^\\?", "", query), "&", fixed = TRUE)[[1]])
  !any(keys %in% c("_state_id_", "_inputs_", "_values_", "showcase"))
}

cache_static_app_ui <- function(app) {
  handler <- app$httpHandler
  cached <- NULL
  app$httpHandler <- function(req) {
    # The UI is a static tag list, so its HTML does not depend on the query string (?go=FALSE, ?dff=..&gene=.. are
    # read in the browser); only Shiny bookmarking makes the page URL-dependent.
    plain <- identical(req$REQUEST_METHOD, "GET") && identical(req$PATH_INFO, "/") &&
      rc_ui_query_cacheable(req$QUERY_STRING %||% "")
    if (!plain) return(handler(req))
    if (!is.null(cached)) return(cached)
    response <- handler(req)
    if (!is.null(response) && isTRUE(response$status == 200L)) cached <<- response
    response
  }
  app
}
