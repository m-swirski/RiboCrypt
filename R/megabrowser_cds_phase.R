# Phase at each displayed base, using the selected transcript's spliced CDS.
# Outside the CDS, the displayed-coordinate phase is extended from the first
# coding overlap; only coding bases have a defined biological CDS frame.
megabrowser_cds_phase <- function(display, cds) {
  if (length(display) != 1L || length(cds) != 1L) return(NULL)
  d <- display[[1L]]; c <- cds[[1L]]
  if (!length(d) || !length(c)) return(NULL)
  ds <- as.character(GenomicRanges::strand(d)); cs <- as.character(GenomicRanges::strand(c))
  if (length(unique(c(ds, cs))) != 1L || !ds[1L] %in% c('+', '-')) return(NULL)
  hits <- GenomicRanges::findOverlaps(d, c, ignore.strand = FALSE)
  if (!length(hits)) return(NULL)
  doff <- cumsum(IRanges::width(d)) - IRanges::width(d)
  coff <- cumsum(IRanges::width(c)) - IRanges::width(c)
  n <- sum(IRanges::width(d)); phase <- rep(NA_integer_, n)
  anchor <- NULL
  for (i in seq_along(hits)) {
    di <- S4Vectors::queryHits(hits)[i]; ci <- S4Vectors::subjectHits(hits)[i]
    lo <- max(IRanges::start(d)[di], IRanges::start(c)[ci])
    hi <- min(IRanges::end(d)[di], IRanges::end(c)[ci])
    genomic <- if (ds[1L] == '+') seq.int(lo, hi) else seq.int(hi, lo)
    dp <- doff[di] + if (ds[1L] == '+') genomic - IRanges::start(d)[di] + 1L else IRanges::end(d)[di] - genomic + 1L
    cp <- coff[ci] + if (ds[1L] == '+') genomic - IRanges::start(c)[ci] + 1L else IRanges::end(c)[ci] - genomic + 1L
    p <- (cp - 1L) %% 3L
    if (any(!is.na(phase[dp]) & phase[dp] != p)) return(NULL)
    phase[dp] <- p
    if (is.null(anchor)) anchor <- (cp[1L] - dp[1L]) %% 3L
  }
  missing <- is.na(phase)
  phase[missing] <- (which(missing) - 1L + anchor) %% 3L
  as.integer(phase)
}
