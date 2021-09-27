multiOmicsController <- function(){
  with(caller_env(),{
seqlevels(target_range) <- seqlevels(annotation)
target_range <- GRangesList(target_range)

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

if (length(proportions) == 0) proportions <- 1
if (!(length(proportions)  %in% c(1, length(reads)))) stop("length of proportions must be 0, 1 or the same as reads list")
if (length(proportions) == 1) proportions <- rep(proportions, length(reads))

if (!display_sequence) {
  max_sum <- 1  - sum(0.03,0.1)
  proportions <- proportions / (sum(proportions) / max_sum)
  proportions <- c(proportions, c(0.03,0.1))
} else {
  max_sum <- 1  - sum(0.035, 0.03,0.1)
  proportions <- proportions / (sum(proportions) / max_sum)
  proportions <- c(proportions, c(0.035, 0.03,0.1))
}
}
)
}
