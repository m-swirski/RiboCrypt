getRiboProfile <- function(grl, footprints, kmer = 1) {
  if (kmer == 1) {
    footprints <- coveragePerTiling(grl, subsetByOverlaps(footprints, grl), as.data.table = TRUE, withFrames=TRUE)
    footprints$frame <- as.factor(footprints$frame)
  } else{
    extended_range <- grl %>% extendLeaders(kmer * 3) %>% extendTrailers(kmer * 3)
    footprints <- coveragePerTiling(extended_range, subsetByOverlaps(footprints, extended_range), as.data.table = TRUE, withFrames=TRUE)
    footprints$frame <- as.factor(footprints$frame)
    footprints$count <- as.numeric(footprints$count)
    footprints$genes <- as.factor(footprints$genes)
    footprints <- footprints[,count := frollmean(count, kmer, fill = 0, align = "center"), by = list(genes,frame)]
    footprints <- footprints[(kmer*3 + 1):(nrow(footprints) - kmer*3)]
    footprints$position <- 1:nrow(footprints)
  }
  return(footprints)
}
