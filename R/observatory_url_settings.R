#' Controls shared by Observatory URL capture, restore and autostart readiness.
#' @noRd
observatory_url_controls <- function() {
  list(
    select = c("gene", "tx", "frames_type", "colors", "summary_track_type", "frames_subset"),
    slider = "kmer",
    numeric = c("extendLeaders", "extendTrailers", "collapsed_introns_width"),
    switch = c("viewMode", "other_tx", "collapsed_introns"),
    checkbox = c("add_uorfs", "add_translon", "add_translons_transcode", "log_scale",
                 "log_scale_protein", "phyloP", "mapability", "withFrames", "summary_track"),
    text = c("genomic_region", "zoom_range", "customSequence", "y_range")
  )
}

#' @noRd
observatory_capture_browser_settings <- function(input) {
  fields <- unlist(observatory_url_controls(), use.names = FALSE)
  settings <- stats::setNames(lapply(fields, function(field) input[[field]]), fields)
  settings$frames_subset <- input$frames_subset %||% character()
  c(settings[!vapply(settings, is.null, logical(1))], list(go = TRUE))
}

#' Keep restoration and autostart aligned with validated annotation defaults.
#' @noRd
observatory_resolve_browser_url_state <- function(state, browser_options) {
  if (is.null(state) || !identical(state$view, "browser")) return(state)
  for (field in c("gene", "tx")) {
    option <- if (field == "gene") "default_gene_meta" else "default_isoform_meta"
    value <- unname(browser_options[option])
    if (shiny::isTruthy(value)) state$browser[[field]] <- value
  }
  state
}

#' Load URL collection annotations before resolving gene and transcript defaults.
#' @noRd
observatory_url_collection_init <- function(state, experiments, initial, exps_dir) {
  experiment <- state$exp
  if (!shiny::isTruthy(experiment) || !experiment %in% experiments ||
      identical(experiment, name(initial))) return(NULL)
  df <- get_exp(experiment, experiments, .GlobalEnv, exps_dir)
  list(df = df, names = get_gene_name_categories(df),
       tx = loadRegion(df), cds = loadRegion(df, "cds"))
}

#' Normalize JSON controls, rejecting malformed or impossible scalar values.
#' @noRd
observatory_normalize_browser_settings <- function(settings) {
  if (is.null(settings)) return(list())
  if (!is.list(settings)) stop("Browser settings must be an object")
  controls <- observatory_url_controls()
  types <- rep(names(controls), lengths(controls))
  names(types) <- unlist(controls, use.names = FALSE)
  types <- c(types, go = "switch")
  fields <- intersect(names(settings), names(types))
  stats::setNames(lapply(fields, function(field) {
    observatory_normalize_control(settings[[field]], field, types[[field]])
  }), fields)
}

#' @noRd
observatory_normalize_control <- function(value, field, type) {
  if (is.null(value)) return(NULL)
  value <- unlist(value, use.names = FALSE)
  if (field == "frames_subset") {
    if (!all(value %in% c("all", "red", "green", "blue"))) stop("Invalid frame subset")
    return(as.character(value))
  }
  if (length(value) != 1L || is.na(value)) stop("Invalid URL control: ", field)
  if (type %in% c("switch", "checkbox")) value <- as.logical(value)
  else if (type %in% c("numeric", "slider")) value <- as.numeric(value)
  else value <- as.character(value)
  if (is.na(value) || (is.numeric(value) && !is.finite(value))) stop("Invalid URL control: ", field)
  if (field == "kmer" && (value < 1 || value > 20)) stop("Invalid kmer")
  if (field %in% c("frames_type", "summary_track_type") &&
      !value %in% c("lines", "columns", "stacks", "area", "heatmap", "animate")) stop("Invalid display type")
  if (field == "colors" && !value %in% c("R", "Color_blind")) stop("Invalid color theme")
  if (field == "y_range") browser_parse_y_range(value)
  value
}

#' Validate the versioned URL object before it reaches reactive UI code.
#' @noRd
observatory_validate_url_state <- function(state) {
  if (!is.list(state) || is.null(names(state))) stop("Invalid URL state")
  version <- state[["v"]]
  if (!is.null(version) && (length(version) != 1L || is.na(version) || version != 1)) stop("Unsupported URL version")
  for (field in c("exp", "view")) {
    value <- state[[field]]
    if (!is.null(value) && (!is.character(value) || length(value) != 1L || is.na(value))) stop("Invalid ", field)
  }
  state$browser <- observatory_normalize_browser_settings(state$browser)
  if (!is.null(state$selections)) observatory_validate_url_selections(state$selections)
  colors <- unlist(state$color_by, use.names = FALSE)
  if (!is.null(colors) && !is.character(colors)) stop("Invalid color columns")
  state
}

#' @noRd
observatory_validate_url_selections <- function(selections) {
  if (!is.list(selections)) stop("Invalid selections")
  ids <- unlist(selections$order, use.names = FALSE)
  if (!is.character(ids) || !length(ids) || anyNA(ids) ||
      any(!grepl("^[1-9][0-9]*$", ids)) || anyDuplicated(ids)) stop("Invalid selection IDs")
  active <- selections$active
  if (length(active) != 1L || !active %in% ids) stop("Invalid active selection")
  runs <- unlist(selections$runs, use.names = FALSE)
  if (!is.null(runs) && (!is.character(runs) || anyNA(runs))) stop("Invalid runs")
  for (field in c("p", "d", "labels")) {
    values <- selections[[field]]
    if (!is.null(values) && !is.list(values)) stop("Invalid selection map")
    for (id in names(values)) {
      value <- unlist(values[[id]], use.names = FALSE)
      if (field == "labels") {
        if (!is.character(value) || length(value) != 1L || is.na(value)) stop("Invalid label")
      } else if (length(value) && (!is.numeric(value) || anyNA(value) ||
                  any(!is.finite(value) | value != trunc(value) | value < 1 | value > length(runs)))) {
        stop("Invalid run indices")
      }
    }
  }
  invisible(NULL)
}

#' Autostart must use restored cohorts, not a nonempty default All merged group.
#' @noRd
observatory_url_selections_ready <- function(expected, actual) {
  if (is.null(expected)) return(TRUE)
  ids <- expected$index
  identical(as.character(names(actual)), as.character(ids)) &&
    all(vapply(ids, function(id) {
      setequal(actual[[id]], expected$data_table_selections[[id]])
    }, logical(1)))
}

#' @noRd
observatory_browser_settings_ready <- function(settings, input) {
  fields <- setdiff(names(settings), "go")
  all(vapply(fields, function(field) {
    expected <- settings[[field]]
    if (is.null(expected)) return(TRUE)
    actual <- input[[field]]
    if (field == "frames_subset") return(setequal(actual, expected))
    length(actual) == 1L && identical(as.character(actual), as.character(expected))
  }, logical(1)))
}
