getRiboProfile <- function(grl, footprints, kmers = 1, kmers_type = "mean") {
  count <- NULL; genes <- NULL # Avoid data.table warning
  not_coverage <- is(footprints, "GenomicRanges") |
    is(footprints, "GAlignments") | is(footprints, "GAlignmentPairs")
  if (kmers == 1) {
    footprints <- coveragePerTiling(grl,
                                    if(not_coverage){
                                      subsetByOverlaps(footprints, grl)}
                                    else {footprints},
                                    as.data.table = TRUE,
                                    withFrames=TRUE, is.sorted = TRUE)
    footprints$frame <- as.factor(footprints$frame)
  } else{
    extended_range <- grl %>% extendLeaders(kmers * 3) %>% extendTrailers(kmers * 3)
    footprints <- coveragePerTiling(extended_range,
                                    if(not_coverage) {
                                      subsetByOverlaps(footprints, grl)}
                                    else {footprints},
                                    as.data.table = TRUE, withFrames=TRUE, is.sorted = TRUE)
    footprints[, frame := as.factor(frame)]
    footprints[, genes := as.factor(genes)]

    footprints <- smoothenSingleSampCoverage(footprints, kmers)
  }
  return(footprints)
}

smoothenSingleSampCoverage <- function(dt, kmer, kmers_type = "mean") {
  dt[, count := as.numeric(count)]
  dt[,count := get(paste("froll", kmers_type, sep = ""))(count, kmer, fill = 0, align = "center"), by = list(genes, frame)]
  dt <- dt[(kmer*3 + 1):(nrow(dt) - kmer*3)]
  dt[, position := 1:nrow(dt)]
  return(dt)
}

smoothenMultiSampCoverage <- function(dt, kmer, kmers_type = "mean", split_by_frame = FALSE) {
  dt[, count := as.numeric(count)]
  if (split_by_frame) {
    if (!("frame" %in% colnames(dt))) {
      dt[, frame := rep(seq(0, 2), length.out = nrow(dt))]
    }
    return(dt[,count := get(paste("froll", kmers_type, sep = ""))(count, kmer, fill = 0, align = "center"), by = .(library, frame)])
  }
  return(dt[,count := get(paste("froll",kmers_type, sep = ""))(count, kmer, fill = 0, align = "center"), by = library])
}

#' Get coverage profile
#'
#' @param grl a GRangesList
#' @param reads GRanges
#' @param kmers 1
#' @param kmers_type "mean"
#' @return data.table of coverage
#' @importFrom dplyr `%>%`
#' @keywords internal
getCoverageProfile <- function(grl, reads, kmers = 1, kmers_type = "mean") {
  not_coverage <- is(reads, "GenomicRanges") |
    is(reads, "GAlignments") | is(reads, "GAlignmentPairs")
  if (kmers == 1) {
    coverage <- coveragePerTiling(grl,
                                  if(not_coverage) {
                                    subsetByOverlaps(reads, grl)}
                                  else {reads},
      as.data.table = TRUE, is.sorted = TRUE)
  } else {
    count <- NULL; genes <- NULL # Avoid data.table warning
    extended_range <- grl %>% extendLeaders(kmers * 3) %>% extendTrailers(kmers * 3)
    coverage <- coveragePerTiling(extended_range,
                                  if(not_coverage) {
                                    subsetByOverlaps(reads, grl)}
                                  else {reads},
                                  as.data.table = TRUE, is.sorted = TRUE)
    coverage$count <- as.numeric(coverage$count)
    coverage$genes <- as.factor(coverage$genes)
    coverage <- coverage[,count := get(paste("froll",kmers_type, sep = ""))(count, kmers, fill = 0, align = "center"), by = genes]
    coverage <- coverage[(kmers*3 + 1):(nrow(coverage) - kmers*3)]
    coverage[, position := seq.int(nrow(coverage))]
  }
  return(coverage)
}

getProfileWrapper <- function(display_range, reads, withFrames, kmers = 1, log_scale = FALSE,
                              kmers_type = "mean", type = "lines", frames_subset = "all") {

  if (withFrames) {
    profile <- getRiboProfile(display_range, reads, kmers, kmers_type = kmers_type)
    if (frames_subset != "all") {
      color_options <- c("red", "green", "blue")
      frame_to_use <- which(color_options %in% frames_subset) - 1
      profile <- profile[frame != frame_to_use, count := 0]
    }
  } else {
    profile <- getCoverageProfile(display_range, reads, kmers, kmers_type = kmers_type)
  }
  if (log_scale) profile[, count := log2(count + 1)]
  return(profile)
}


getProfileAnimate <- function(display_range, reads, withFrames, kmers = 1, kmers_type = "mean") {
  if (withFrames) {
    profile <- getRiboProfile(display_range, reads, kmers, kmers_type = kmers_type)

  } else {
    not_coverage <- is(reads, "GenomicRanges") |
      is(reads, "GAlignments") | is(reads, "GAlignmentPairs")
    profile <- coveragePerTiling(display_range, if(not_coverage) {
      subsetByOverlaps(reads, display_range)}
      else {reads}, as.data.table = TRUE)
  }
  profile
}
