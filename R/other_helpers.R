
#' Get index
#' @param ref_granges a GRanges object
#' @return integer vector, indices
#' @importFrom BiocGenerics start end
#' @keywords internal
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

#' Match multiple patterns
#' @param patterns character
#' @param Seq a DNAStringSet
#' @return integer vector, indices (named with pattern hit)
#' @importFrom Biostrings matchPattern
#' @keywords internal
matchMultiplePatterns <- function(patterns, Seq) {
  matches <- c()
  for (pattern in patterns) {
    cur_matches <- matchPattern(pattern, Seq, fixed = FALSE)
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

#' Match to GRanges
#' @param matches integer vector, indices
#' @param ref_granges GRanges
#' @return GRanges object
#' @importFrom GenomeInfoDb seqnames
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @keywords internal
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

#' Get antisense
#'
#' @importFrom BiocGenerics strand strand<-
#' @return a GRangesList
#' @keywords internal
antisense <- function(grl) {
  gr <- unlistGrl(grl)
  pos_pointer <- strandBool(gr)
  strand(gr[pos_pointer]) <- "-"
  strand(gr[!pos_pointer]) <- "+"
  grl <- groupGRangesBy(gr)
  return(grl)
}

#' Trim overlaps
#'
#' @param overlaps GRanges
#' @param display_range GRanges
#' @return GRanges
#' @importFrom BiocGenerics start<- end<-
#' @keywords internal
trimOverlaps <- function(overlaps, display_range) {
  flanks <- unlistGrl(display_range)

  start_out_of_bounds <- start(overlaps) < min(start(flanks))
  end_out_of_bounds <- end(overlaps) > max(end(flanks))
  if (any(start_out_of_bounds)) {
    start(overlaps)[start_out_of_bounds] <- min(start(flanks))
  }
  if (any(end_out_of_bounds)) {
    end(overlaps)[end_out_of_bounds] <- max(end(flanks))
  }
  return(overlaps)
}

# selectCols <- function(cols, locations) {
#   matches <- match(names(locations), names(cols))
#   duplications<- duplicated(matches)
#   additions <- cumsum(duplications)
#   additions[!duplications] <- 0
#   matches <- matches + additions
#
#   return(cols[matches])
# }

selectFrames <- function(frames, locations) {
  matches <- match(names(locations), names(frames))
  duplications<- duplicated(matches)
  additions <- cumsum(duplications)
  additions[!duplications] <- 0
  matches <- matches + additions

  return(frames[matches])
}


colour_bars <- function(locations, overlaps, display_range) {
  all_colors <- rep("gray", length(locations))
  names(all_colors) <- names(locations)
  is_cds <- locations$type == "cds"
  if (all(!is_cds)) return(all_colors)
  locations <- locations[is_cds]
  overlaps <- overlaps[overlaps$type == "cds"]
  frames <- start(locations) - 1
  frames <- frames + locations$rel_frame_orf
  colors <- rc_rgb()[frames %% 3 + 1]
  tx_mode <- 1 %in% start(locations)
  if (tx_mode) {
    plus_stranded <- as.logical(strand(display_range[[1]][1]) == "+")
    if (plus_stranded) {
      ov_starts <- start(overlaps)
      target_start <- min(start(display_range))
      frames <- ov_starts - target_start
      frames <- frames + overlaps$rel_frame_exon
      colors_general <- rc_rgb()[frames %% 3 + 1]
    } else {
      ov_starts <- end(overlaps)
      target_start <- max(end(display_range))
      frames <- ov_starts - target_start
      frames <- frames - overlaps$rel_frame_exon
      frames <- - frames
      colors_general <- rc_rgb()[c(3,1,2)][(frames+1) %% 3 + 1]
    }
    is_tx_init <- start(locations) == 1
    colors[is_tx_init] <- colors_general[which(is_tx_init)]
  }
  all_colors[is_cds] <- colors
  return(all_colors)
}

#'
#' @import data.table
#' @importFrom Biostrings width
#' @keywords internal
getRelativeFrames <- function(overlaps) {
  dt <- data.table(names = names(overlaps),
                   width = width(overlaps),
                   type = overlaps$type)
  dt[,cum_width := cumsum(width), .(names, type)]
  dt[, rel_frame := c(0,-cum_width %% 3)[1:length(width)], .(names, type)]
  return(dt$rel_frame)
}


realA <- function(y,z,Ea,Eb,Ec) {
  eq <- matrix(c(1 - y - z,z,y,
                 y,1 - y - z, z,
                 z,y,1 - y - z), 3,3, byrow = TRUE
  )
  sol <- c(Ea, Eb, Ec)
  solve(eq,sol)
}

