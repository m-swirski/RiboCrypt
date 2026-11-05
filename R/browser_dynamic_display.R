# Select a rendering strategy from the visible span, in nucleotide coordinates.
#' @noRd
browser_dynamic_style <- function(span) {
  if (length(span) != 1L || !is.finite(span) || span <= 0) return(list(type = "columns", kmer = 1L))
  if (span > 10000) return(list(type = "lines", kmer = 9L))
  if (span > 500) return(list(type = "area", kmer = 9L))
  list(type = "columns", kmer = 1L)
}

#' @noRd
browser_dynamic_context <- function(input) {
  paste(vapply(c("go", "tx", "genomic_region", "zoom_range", "viewMode", "extendLeaders", "extendTrailers", "collapsed_introns", "collapsed_introns_width"),
    function(nm) paste(input[[nm]] %||% "", collapse = ":"), character(1)), collapse = "|:|")
}

#' @noRd
browser_dynamic_settings <- function(input, display_region, zoom_range) {
  active <- identical(input$display_mode %||% "manual", "dynamic")
  if (!active) return(list(active = FALSE, type = input$frames_type, kmer = input$kmer, range = NULL, context = NULL))
  context <- browser_dynamic_context(input)
  extension <- function(value) {
    value <- suppressWarnings(as.numeric(value))
    if (length(value) != 1L || !is.finite(value)) 0 else value
  }
  display_length <- as.numeric(widthPerGroup(display_region, FALSE))[1L] +
    extension(input$extendLeaders) + extension(input$extendTrailers)
  if (!is.finite(display_length) || display_length <= 0) display_length <- 1
  range <- if (length(zoom_range) == 2L && all(is.finite(zoom_range))) zoom_range else c(1, display_length)
  event <- input$dynamic_view
  if (is.list(event) && identical(event$context, context)) {
    candidate <- suppressWarnings(as.numeric(unlist(event$range)))
    if (length(candidate) == 2L && all(is.finite(candidate)) && candidate[2L] > candidate[1L]) range <- candidate
  }
  style <- browser_dynamic_style(diff(range) + 1)
  list(active = TRUE, type = style$type, kmer = style$kmer, range = range, context = context)
}
