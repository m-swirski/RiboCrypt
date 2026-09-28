#' Parse an automatic range, a maximum, or a minimum:maximum pair.
#' @noRd
browser_parse_y_range <- function(value = "auto") {
  if (is.null(value)) return(NULL)
  if (length(value) != 1L || is.na(value)) stop("Y-axis range must be a single value.", call. = FALSE)
  value <- tolower(trimws(as.character(value)))
  if (identical(value, "auto")) return(NULL)
  number <- "(?:[0-9]+(?:\\.[0-9]*)?|\\.[0-9]+)(?:e[+-]?[0-9]+)?[km]?"
  if (!grepl(paste0("^", number, "(?:\\s*:\\s*", number, ")?$"), value, perl = TRUE)) {
    stop("Y-axis range: enter 'auto', a positive maximum (500000 or 500k), or min:max (0:500000).", call. = FALSE)
  }
  limits <- vapply(strsplit(value, ":", fixed = TRUE)[[1]], browser_y_range_number, numeric(1))
  if (length(limits) == 1L) limits <- c(0, limits)
  if (any(!is.finite(limits)) || limits[1] >= limits[2]) {
    stop("Y-axis limits must be finite, non-negative numbers, with maximum greater than minimum.", call. = FALSE)
  }
  unname(limits)
}

#' @noRd
browser_y_range_number <- function(value) {
  value <- trimws(value)
  multiplier <- if (endsWith(value, "k")) 1e3 else if (endsWith(value, "m")) 1e6 else 1
  suppressWarnings(as.numeric(sub("[km]$", "", value))) * multiplier
}

#' Fail early, before annotation or coverage loading, with a user-facing message.
#' @noRd
browser_validate_y_range <- function(value) {
  tryCatch(browser_parse_y_range(value), error = function(error) {
    message <- conditionMessage(error)
    shiny::showModal(shiny::modalDialog(title = "Invalid Y-axis range", message, easyClose = TRUE))
    shiny::validate(shiny::need(FALSE, message))
  })
}

#' Override only coverage axes; annotation and heatmap y coordinates are not counts.
#' @noRd
browser_apply_y_range <- function(plot, limits, track_types) {
  if (is.null(limits)) return(plot)
  axes <- which(track_types != "heatmap")
  for (index in axes) {
    name <- if (index == 1L) "yaxis" else paste0("yaxis", index)
    patch <- list(range = limits, autorange = FALSE, tickmode = "auto", nticks = 5,
                  tickvals = NULL, ticktext = NULL)
    plot$x$layout[[name]] <- utils::modifyList(plot$x$layout[[name]], patch, keep.null = TRUE)
    plot <- do.call(plotly::layout, c(list(plot), stats::setNames(list(patch), name)))
  }
  plot
}

#' Match subplot ordering, including the optional summary and animation panel.
#' @noRd
browser_coverage_track_types <- function(controls, profiles) {
  count <- if (identical(controls$frames_type, "animate")) 1L else length(profiles)
  types <- rep(controls$frames_type, count)
  if (isTRUE(controls$summary_track)) types <- c(controls$summary_track_type, rev(types))
  types
}
