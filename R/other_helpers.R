getIndexes <- function(ref_granges) {
  if (nrun(strand(ref_granges)) > 1) stop("More than one strand!")
  if (nrun(seqnames(ref_granges)) > 1) stop("More than one seqname!")
  starts <- start(ref_granges)
  ends <- end(ref_granges)
  if (as.vector(strand(ref_granges))[1] == "+"){
    indexes <- mapply(function(x,y) x:y , starts, ends) %>% unlist %>% sort
  }else {
    indexes <- mapply(function(x,y) x:y , starts, ends) %>% unlist %>% sort(decreasing = TRUE)
  }
  return(indexes)
}
matchMultiplePatterns <- function(patterns, Seq) {
  matches <- c()
  for (pattern in patterns) {
    cur_matches <- matchPattern(pattern, Seq)
    cur_matches <- cur_matches@ranges@start
    names(cur_matches) <- rep(pattern, length(cur_matches))
    matches <- c(matches, cur_matches)
  }
  return(matches)
}

findPatterns <- function(patterns, sequence, frame = c(0,1,2), max_dist = NULL, min_dist = NULL) {
  matches <- matchMultiplePatterns(patterns = patterns, Seq = sequence)
  if (!identical(frame, c(0,1,2))){
    match_frames <- (matches - 1) %% 3
    match_frames <- match_frames %in% frame
    matches <- matches[match_frames]
  }
  if (!is.null(max_dist)){
    matches <- matches[matches <= max_dist]
  }
  if (!is.null(min_dist)) {
    matches <- matches[matches >= min_dist]
  }
  return(matches)
}

matchToGRanges <- function(matches, ref_granges) {
  indexes <- getIndexes(ref_granges = ref_granges)
  starts <- indexes[matches]
  match_ranges <- IRanges(start = starts, end = starts)
  gene_strand = strand(ref_granges)[1]
  gene_seqname = seqnames(ref_granges)[1]
  match_granges <- GRanges(seqnames = gene_seqname, strand = gene_strand, ranges = match_ranges)
  return(match_granges)
}

matchWrapper <- function(patterns, ref_granges,sequence, frame = c(0,1,2), max_dist = NULL, min_dist = NULL){
  matches <- findPatterns(patterns = patterns, sequence = sequence, frame = frame, max_dist = max_dist, min_dist = min_dist)
  if (length(matches) != 0) {
    match_granges <- matchToGRanges(matches = matches, ref_granges = ref_granges)
    return(match_granges)
  }
}


antisense <- function(grl) {
  gr <- unlistGrl(grl)
  pos_pointer <- strandBool(gr)
  strand(gr[pos_pointer]) <- "-"
  strand(gr[!pos_pointer]) <- "+"
  grl <- groupGRangesBy(gr)
  return(grl)

}


trimOverlaps <- function(overlaps, target_range) {
  tr <- unlistGrl(target_range)
  start_indices <- start(overlaps) < min(start(tr))
  end_indices <- end(overlaps) > max(end(tr))
  if (TRUE %in% start_indices) {
    start(overlaps)[start(overlaps) < min(start(tr))] <- min(start(tr))

  }
  if (TRUE %in% end_indices) {
    end(overlaps)[end(overlaps) > max(end(tr))] <- max(end(tr))
  }
  return(overlaps)
}

selectCols <- function(cols, locations) {
  matches <- match(names(locations), names(cols))
  duplications<- duplicated(matches)
  additions <- cumsum(duplications)
  additions[!duplications] <- 0
  matches <- matches + additions

  return(cols[matches])
}


colour_bars <- function(overlaps, target_range) {
  if (as.logical(strand(target_range[[1]][1]) == "+")) {
    ov_starts <- start(overlaps)
    target_start <- min(start(target_range))
    frames <- ov_starts - target_start
    frames <- frames + overlaps$rel_frame
    colors <- c("#F8766D","#00BA38","#619CFF")[frames %% 3 + 1]
    return(colors)
  } else{
    ov_starts <- end(overlaps)
    target_start <- max(end(target_range))
    frames <- ov_starts - target_start
    frames <- frames - overlaps$rel_frame
    frames <- - frames
    colors <- c("#619CFF","#F8766D","#00BA38")[(frames+1) %% 3 + 1]
    return(colors)
  }
}

#'
#' @import data.table
#' @keywords internal
getRelativeFrames <- function(overlaps) {
  dt <- data.table(names = names(overlaps),
                   width = width(overlaps))
  dt[,cum_width := cumsum(width), names]
  dt[, rel_frame := c(0,-cum_width %% 3)[1:length(width)], names]
  return(dt$rel_frame)
}
