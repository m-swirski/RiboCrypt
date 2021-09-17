getRiboProfile <- function(grl, footprints, kmer = 1) {
  if (kmer == 1) {
    footprints <- coveragePerTiling(grl, subsetByOverlaps(footprints, grl), as.data.table = TRUE, withFrames=TRUE, is.sorted = T)
    footprints$frame <- as.factor(footprints$frame)
  } else{
    extended_range <- grl %>% extendLeaders(kmer * 3) %>% extendTrailers(kmer * 3)
    footprints <- coveragePerTiling(extended_range, subsetByOverlaps(footprints, extended_range), as.data.table = TRUE, withFrames=TRUE, is.sorted = T)
    footprints$frame <- as.factor(footprints$frame)
    footprints$count <- as.numeric(footprints$count)
    footprints$genes <- as.factor(footprints$genes)
    footprints <- footprints[,count := frollmean(count, kmer, fill = 0, align = "center"), by = list(genes,frame)]
    footprints <- footprints[(kmer*3 + 1):(nrow(footprints) - kmer*3)]
    footprints$position <- 1:nrow(footprints)
  }
  return(footprints)
}
getCoverageProfile <- function(grl, reads, kmer = 1) {
  if (kmer == 1) {
    coverage <- coveragePerTiling(grl, subsetByOverlaps(reads, grl), as.data.table = TRUE, is.sorted = T)
  } else {
    extended_range <- grl %>% extendLeaders(kmer * 3) %>% extendTrailers(kmer * 3)
    coverage <- coveragePerTiling(extended_range, subsetByOverlaps(reads, extended_range), as.data.table = TRUE, is.sorted = T)
    coverage$count <- as.numeric(coverage$count)
    coverage$genes <- as.factor(coverage$genes)
    coverage <- coverage[,count := frollmean(count, kmer, fill = 0, align = "center"), by = genes]
    coverage <- coverage[(kmer*3 + 1):(nrow(coverage) - kmer*3)]
    coverage$position <- 1:nrow(coverage)
  }
  return(coverage)
}

getStackProfile <- function(grl, footprints, kmer) {
  profile <- getRiboProfile(grl, footprints, kmer)
  maxpos <- max(profile$position)
  profile <- profile[,.(count = rep(count, 3), new_position = c(position, position + 1, position + 2)),.(position,frame)]
  profile$position <- profile$new_position
  profile <- profile[position <= maxpos]
  return(profile)
}


getProfileAnimate <- function(target_range, reads, withFrames, kmer = 1) {
  if (withFrames) {
    profile <- getRiboProfile(target_range, reads, kmer)

  } else {
    profile <- coveragePerTiling(target_range, subsetByOverlaps(reads, target_range), as.data.table = TRUE)
  }
  profile
}
