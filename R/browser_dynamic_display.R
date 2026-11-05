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

# Compact sidecar: x coordinates and raw counts already exist in column traces.
#' @noRd
browser_dynamic_client_config <- function(raw, smooth, frame_colors, input_id, status_id) {
  tracks <- lapply(seq_along(raw), function(i) {
    levels <- track_frame_levels(raw[[i]])
    colors <- track_frame_colors(levels, frame_colors)
    fills <- track_frame_colors(levels, frame_colors, with_alpha = TRUE)
    frames <- lapply(levels, function(level) {
      r <- raw[[i]][as.character(frame) == level]
      s <- smooth[[i]][as.character(frame) == level]
      if (!identical(r$position, s$position)) stop("Dynamic coverage position mismatch")
      list(y = unname(s$count), color = unname(colors[[level]]),
           fill = grDevices::adjustcolor(fills[[level]], alpha.f = track_area_fill_alpha()))
    })
    names(frames) <- levels
    frames
  })
  list(tracks = tracks, input_id = input_id, status_id = status_id,
       simplify = isTRUE(getOption("RiboCrypt.browser_simplify", TRUE)))
}

# Load the padded RFP coverage once, then derive the exact k=1 and k=9
# versions. The 27 nt padding is a multiple of three, preserving frame phase.
# Fall back to the existing readers if the complete padding is unavailable.
#' @noRd
browser_dynamic_dual_profiles <- function(display_range, reads, frames_subset, log_scale) {
  extended <- extendTrailers(extendLeaders(display_range, 27L), 27L)
  original_length <- as.numeric(widthPerGroup(display_range, FALSE))[1L]
  padded_length <- as.numeric(widthPerGroup(extended, FALSE))[1L]
  if (!is.finite(original_length) || !is.finite(padded_length) || padded_length != original_length + 54L) return(NULL)
  per_track <- lapply(reads, function(read) {
    padded <- getRiboProfile(extended, read, kmers = 1L)
    if (nrow(padded) != padded_length) return(NULL)
    raw <- data.table::copy(padded[seq.int(28L, 27L + original_length)])
    raw[, position := seq_len(.N)]
    smooth <- data.table::copy(padded)
    smooth[, genes := as.factor(genes)]
    smooth <- smoothenSingleSampCoverage(smooth, 9L)
    if (all(frames_subset != "all")) {
      use <- which(c("red", "green", "blue") %in% frames_subset) - 1L
      raw[!(frame %in% use), count := 0]
      smooth[!(frame %in% use), count := 0]
    }
    if (isTRUE(log_scale)) {
      raw[, count := log2(count + 1)]
      smooth[, count := log2(count + 1)]
    }
    list(raw = raw, smooth = smooth)
  })
  if (any(vapply(per_track, is.null, logical(1)))) return(NULL)
  list(raw = lapply(per_track, `[[`, "raw"), smooth = lapply(per_track, `[[`, "smooth"))
}

#' @noRd
browser_dynamic_client_supported <- function(active, with_frames, summary_track) {
  isTRUE(active) && length(with_frames) > 0L && !anyNA(with_frames) &&
    all(with_frames) && !isTRUE(summary_track)
}
