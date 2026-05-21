
dt_fft <- function(dt) {
  read_lengths <- unique(dt$fraction)
  fft_dt <- data.table()
  for (i in read_lengths) {
    spec <- spec.pgram(x = dt[fraction == i, ]$score, plot = FALSE)
    fft_dt <- rbindlist(list(fft_dt, data.table(read_length = i, 
                                                amplitude = spec$spec, periods = 1/spec$freq)))
  }
  return(fft_dt)
}

#' Tile coverage while keeping exon identities
#'
#' `coveragePerTiling()` concatenates the ranges in each `GRangesList`
#' element. This helper keeps that efficient coverage call, then annotates each
#' tiled row with the source exon so exon-level transcript QC can be done
#' without a second coverage pass.
#'
#' @param grl a `GRangesList` with one transcript or feature per element and
#'   one or more exon-like ranges inside each element.
#' @param reads reads or coverage input accepted by
#'   [ORFik::coveragePerTiling()].
#' @param ... extra arguments passed to [ORFik::coveragePerTiling()].
#'
#' @return a `data.table` from [ORFik::coveragePerTiling()] with
#'   `transcript_id`, `transcript_index`, `exon_id`, `exon_rank`,
#'   `exon_width`, and `exon_position` columns.
#' @export
coveragePerExonTiling <- function(grl, reads, ...) {
  stopifnot(is(grl, "GRangesList"))

  coverage_call <- function(target) {
    coveragePerTiling(target, reads, as.data.table = TRUE, ...)
  }
  coverage <- tryCatch(
    coverage_call(grl),
    error = function(e) {
      split_bigwig_message <- grepl(
        "bigwig random access only supported for single gene",
        conditionMessage(e),
        fixed = TRUE
      )
      if (!split_bigwig_message || length(grl) <= 1L) stop(e)
      rbindlist(lapply(seq_along(grl), function(i) {
        per_tx <- as.data.table(coverage_call(grl[i]))
        per_tx[, genes := i]
        per_tx
      }), fill = TRUE)
    }
  )
  coverage <- as.data.table(coverage)
  if (nrow(coverage) == 0) {
    coverage[, `:=`(
      transcript_id = character(),
      transcript_index = integer(),
      exon_id = character(),
      exon_rank = integer(),
      exon_width = integer(),
      exon_position = integer()
    )]
    return(coverage[])
  }

  if (!all(c("genes", "position") %in% names(coverage))) {
    stop("coveragePerTiling() must return genes and position columns.")
  }

  grl_names <- names(grl)
  if (is.null(grl_names)) grl_names <- rep("", length(grl))
  unnamed <- is.na(grl_names) | !nzchar(grl_names)
  grl_names[unnamed] <- paste0("group_", which(unnamed))

  exon_map <- rbindlist(lapply(seq_along(grl), function(i) {
    exon_widths <- as.integer(width(grl[[i]]))
    if (length(exon_widths) == 0) return(NULL)
    exon_ends <- cumsum(exon_widths)
    exon_starts <- c(1L, head(exon_ends, -1L) + 1L)
    data.table(
      transcript_index = i,
      transcript_id = grl_names[[i]],
      exon_rank = seq_along(exon_widths),
      exon_start_position = exon_starts,
      exon_end_position = exon_ends,
      exon_width = exon_widths
    )
  }), fill = TRUE)
  exon_map[, exon_id := paste0(transcript_id, ":exon", exon_rank)]

  coverage[, transcript_index := as.integer(genes)]
  coverage <- merge(
    coverage,
    exon_map,
    by = "transcript_index",
    allow.cartesian = TRUE,
    sort = FALSE
  )
  coverage <- coverage[
    position >= exon_start_position & position <= exon_end_position
  ]
  coverage[, exon_position := position - exon_start_position + 1L]
  coverage[, c("exon_start_position", "exon_end_position") := NULL]
  setcolorder(
    coverage,
    c(
      "count", "genes", "position", "transcript_id", "transcript_index",
      "exon_id", "exon_rank", "exon_width", "exon_position"
    )
  )
  coverage[]
}

#' Flag CDS exons with outlier-low coverage
#'
#' This QC is intended for a CDS `GRangesList` and a broad aggregate read
#' track. A transcript is flagged when a measurable exon has mean coverage
#' below a chosen fraction of the transcript's mean exon coverage.
#'
#' @param grl a CDS `GRangesList` with one transcript per element.
#' @param reads reads or coverage input accepted by
#'   [ORFik::coveragePerTiling()].
#' @param low_exon_fraction coverage fraction below which an exon is marked
#'   low versus the mean exon coverage for its transcript. Default is `0.1`.
#' @param min_transcript_counts minimum tiled transcript counts required before
#'   low-exon flags are emitted. Default is `10`.
#' @param ... extra arguments passed to [coveragePerExonTiling()].
#'
#' @return one `data.table` row per exon with tiled counts, exon mean coverage,
#'   transcript-level reference coverage, per-exon ratio, and transcript flags.
#' @export
cdsExonCoverageQC <- function(grl, reads, low_exon_fraction = 0.1,
                              min_transcript_counts = 10, ...) {
  coverage <- coveragePerExonTiling(grl, reads, ...)
  if (nrow(coverage) == 0) {
    return(data.table(
      transcript_id = character(),
      transcript_index = integer(),
      exon_id = character(),
      exon_rank = integer(),
      exon_width = integer(),
      exon_counts = numeric(),
      exon_mean_coverage = numeric(),
      transcript_exon_mean_coverage = numeric(),
      exon_vs_transcript_mean = numeric(),
      transcript_counts = numeric(),
      n_exons = integer(),
      exon_qc_measurable = logical(),
      low_coverage_exon = logical(),
      transcript_has_low_coverage_exon = logical()
    ))
  }

  count <- NULL
  summary <- coverage[
    ,
    .(
      exon_width = unique(exon_width)[1],
      exon_counts = sum(count, na.rm = TRUE),
      exon_mean_coverage = mean(count, na.rm = TRUE)
    ),
    by = .(transcript_id, transcript_index, exon_id, exon_rank)
  ]
  summary[
    ,
    `:=`(
      transcript_exon_mean_coverage = mean(exon_mean_coverage, na.rm = TRUE),
      transcript_counts = sum(exon_counts, na.rm = TRUE),
      n_exons = .N
    ),
    by = transcript_index
  ]
  summary[, exon_vs_transcript_mean := fifelse(
    is.finite(transcript_exon_mean_coverage) &
      transcript_exon_mean_coverage > 0,
    exon_mean_coverage / transcript_exon_mean_coverage,
    NA_real_
  )]
  summary[, exon_qc_measurable :=
            n_exons > 1L & is.finite(transcript_counts) &
            transcript_counts >= min_transcript_counts &
            is.finite(exon_vs_transcript_mean)]
  summary[, low_coverage_exon :=
            exon_qc_measurable & exon_vs_transcript_mean < low_exon_fraction]
  summary[
    ,
    transcript_has_low_coverage_exon := any(low_coverage_exon, na.rm = TRUE),
    by = transcript_index
  ]
  setorder(summary, transcript_index, exon_rank)
  summary[]
}
