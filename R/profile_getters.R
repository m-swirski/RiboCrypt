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
    footprints$frame <- as.factor(footprints$frame)
    footprints$count <- as.numeric(footprints$count)
    footprints$genes <- as.factor(footprints$genes)
    footprints <- footprints[,count := get(paste("froll",kmers_type, sep = ""))(count, kmers, fill = 0, align = "center"), by = list(genes,frame)]
    footprints <- footprints[(kmers*3 + 1):(nrow(footprints) - kmers*3)]
    footprints$position <- 1:nrow(footprints)
  }
  return(footprints)
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
    coverage$position <- 1:nrow(coverage)
  }
  return(coverage)
}

getStackProfile <- function(grl, footprints, kmers, kmers_type = "mean") {
  count <- NULL # Avoid data.table warning
  profile <- getRiboProfile(grl, footprints, kmers, kmers_type = kmers_type)
  maxpos <- max(profile$position)
  profile <- profile[,.(count = rep(count, 3), new_position = c(position, position + 1, position + 2)),.(position,frame)]
  profile$position <- profile$new_position
  profile <- profile[position <= maxpos]
  return(profile)
}

getProfileWrapper <- function(display_range, reads, withFrames, kmers = 1,
                              kmers_type = "mean", type = "lines") {

  if (withFrames) {
    if (type %in% c("stacks", "area" )) {
      profile <- getStackProfile(display_range, reads, kmers, kmers_type = kmers_type)
    } else profile <- getRiboProfile(display_range, reads, kmers, kmers_type = kmers_type)
  } else {
    profile <- getCoverageProfile(display_range, reads, kmers, kmers_type = kmers_type)
  }
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
