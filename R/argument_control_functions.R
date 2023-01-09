# Controller for input validation
multiOmicsController <- function() {
  with(rlang::caller_env(), {
  # Load ranges if defined as character
  if (is(display_range, "character")) {
    display_range <- loadRegion(df, part = "tx", names.keep = display_range)
    if (length(display_range) == 0) stop("display_range specified as name of non existing transcript,
                                       did you specify a name from a different genome?")
  }
  if (is(annotation, "character")) {
    annotation <- loadRegion(df, part = annotation)
    if (length(annotation) == 0) stop("annotation specified on valid transcript,
                                    but does not contain this region, did you specify cds on
                                    a non coding RNA or leader for mRNA without defined leader?")
  }
  if (viewMode == "genomic") {
    display_range <- flankPerGroup(display_range)
  }
  seqlevels(display_range) <- seqlevels(annotation)
  display_range <- GRangesList(display_range)
  if (leader_extension != 0) display_range <- extendLeaders(display_range, leader_extension)
  if (trailer_extension != 0) display_range <- extendTrailers(display_range, trailer_extension)
  # if (is(annotation, "GRangesList")) annotation <- unlist(annotation)

  if (!is.null(annotation_names)) {
    if (length(annotation_names) == 1){
      if (annotation_names %in% annotation) {
        names(annotation) <- mcols(annotation)[[annotation_names]]
      } else stop(wmsg("wrong annotation_names argument"))
    } else if (length(annotation_names) == length(annotation)) {
      names(annotation) <- annotation_names
    } else stop(wmsg("wrong annotation_names argument"))
  }

  if(!is(reads, "list")) reads <- list(reads)

  if (length(withFrames) == 0) withFrames <- FALSE
  if (!(length(withFrames)  %in% c(1, length(reads)))) stop("length of withFrames must be 0, 1 or the same as reads list")
  if (length(withFrames) == 1) withFrames <- rep(withFrames, length(reads))

  if (length(colors) == 0) colors <- 1:length(reads)
  if (!(length(colors)  %in% c(1, length(reads)))) stop("length of colors must be 0, 1 or the same as reads list")
  if (length(colors) == 1) colors <- rep(colors, length(reads))

  if (length(kmers) == 0) kmers <- 1
  if (!(length(kmers)  %in% c(1, length(reads)))) stop("length of kmers must be 0, 1 or the same as reads list")
  if (length(kmers) == 1) kmers <- rep(kmers, length(reads))

  if (length(ylabels) == 0) ylabels <- as.character(1:length(reads))
  if (!(length(ylabels)  %in% c(1, length(reads)))) stop("length of ylabels must be 0, 1 or the same as reads list")
  if (length(ylabels) == 1) ylabels <- rep(ylabels, length(reads))
  ylabels_full_name <- ylabels
  if (length(ylabels) > 5) ylabels <- 1:length(reads)

  if (length(lib_proportions) == 0) lib_proportions <- 1
  if (!(length(lib_proportions)  %in% c(1, length(reads)))) stop("length of lib_proportions must be 0, 1 or the same as reads list")
  if (length(lib_proportions) == 1) lib_proportions <- rep(lib_proportions, length(reads))
  lib_proportions <- lib_proportions / sum(lib_proportions)
  if (length(annotation_proportions) == 0) {
    if (display_sequence %in% c("none", FALSE)) {
      annotation_proportions <- c(0.35,0.65)
    } else annotation_proportions <- c(0.2,0.2,0.6)
  } else annotation_proportions <- annotation_proportions / sum(annotation_proportions)
  proportions <- c(lib_proportions * lib_to_annotation_proportions[1], annotation_proportions * lib_to_annotation_proportions[2])
  if (summary_track) proportions <- c(0.2, proportions * 0.8)
  }

  )
}
