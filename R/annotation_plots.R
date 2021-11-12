#' Create sequence panel for RiboCrypt
#'
#' @param start_codons character vector, default "ATG"
#' @param stop_codons character vector, default c("TAA", "TAG", "TGA")
#' @param custom_motif character vector, default NULL.
#' @import ggplot2
#' @return a ggplot object
#' @keywords internal
createSeqPanel <- function(sequence, start_codons = "ATG", stop_codons = c("TAA", "TAG", "TGA"), frame = 1,
                           custom_motif = NULL) {
  starts <- matchMultiplePatterns(start_codons, sequence)
  stops <- matchMultiplePatterns(stop_codons, sequence)
  starts_0 <- starts[(starts - 1) %% 3 == 0 ]
  starts_1 <- starts[(starts - 1) %% 3 == 1 ]
  starts_2 <- starts[(starts - 1) %% 3 == 2 ]
  stops_0 <- stops[(stops - 1) %% 3 == 0 ]
  stops_1 <- stops[(stops - 1) %% 3 == 1 ]
  stops_2 <- stops[(stops - 1) %% 3 == 2 ]
  if (!is.null(custom_motif)) {
    custom <- matchMultiplePatterns(custom_motif, sequence)
    custom_0 <- custom[(custom - 1) %% 3 == 0 ]
    custom_1 <- custom[(custom - 1) %% 3 == 1 ]
    custom_2 <- custom[(custom - 1) %% 3 == 2 ]
  }
  suppressWarnings({
  fig <- ggplot() + geom_rect(aes(ymin = 1, ymax = 2, xmin = 1, xmax = length(sequence), frame = frame), fill = "#F8766D") +
    geom_rect(aes(ymin = 0, ymax = 1, xmin = 1, xmax = length(sequence), frame = frame), fill = "#00BA38") +
    geom_rect(aes(ymin = -1, ymax = 0, xmin = 1, xmax = length(sequence), frame = frame), fill = "#619CFF") +
    ylab("frame") +
    xlab("position [nt]") +
    theme(plot.margin = unit(c(0,0,0,0), "pt")) +
    scale_y_continuous(breaks = c(-0.5,0.5, 1.5), labels = c("2","1", "0"), expand = c(0,0)) + scale_x_continuous(expand = c(0,0))
  if (length(starts_0) > 0) fig <- fig + geom_segment(aes(y = 1 , yend = 2, x = starts_0, xend = starts_0), frame = frame, col = "white")
  if (length(stops_0) > 0) fig <- fig + geom_segment(aes(y = 1, yend = 2, x = stops_0, xend = stops_0), frame = frame)
  if (length(starts_1) > 0) fig <- fig + geom_segment(aes(y = 0 , yend = 1, x = starts_1, xend = starts_1), frame = frame, col = "white")
  if (length(stops_1) > 0) fig <- fig + geom_segment(aes(y = 0, yend = 1, x = stops_1, xend = stops_1), frame = frame)
  if (length(starts_2) > 0) fig <- fig + geom_segment(aes(y = -1 , yend = 0, x = starts_2, xend = starts_2), frame=frame, col = "white")
  if (length(stops_2) > 0) fig <- fig + geom_segment(aes(y = -1, yend = 0, x = stops_2, xend = stops_2), frame=frame)
  if (!is.null(custom_motif)) {
    if (length(custom_0) > 0) fig <- fig + geom_segment(aes(y = 1 , yend = 2, x = custom_0, xend = custom_0), frame = frame, col = "purple")
    if (length(custom_1) > 0) fig <- fig + geom_segment(aes(y = 0 , yend = 1, x = custom_1, xend = custom_1), frame = frame, col = "purple")
    if (length(custom_2) > 0) fig <- fig + geom_segment(aes(y = -1 , yend = 0, x = custom_2, xend = custom_2), frame=frame, col = "purple")
  }
  })

  return(fig)
}


#'
#' @keywords internal
overlaps_to_layers <- function(overlaps) {

  v <- unique(overlaps[,2][[1]])

  indexes <- vector(mode = "list", length = length(v))
  names(indexes) <- as.character(v)

  for (x in seq(nrow(overlaps))) {
    indexes[[as.character(overlaps[x,2])]] <-  c(indexes[[as.character(overlaps[x,2])]], get_current_index(indexes[[as.character(overlaps[x,1])]])  )
  }
  sapply(indexes, get_current_index)
}

#'
#' @keywords internal
get_current_index <- function(v) {
  min(which(!(1:(max(c(0,v)+1)) %in% c(0,v))))
}

#'
#' @keywords internal
geneTrackLayer <- function(grl) {

  grl_flanks <- flankPerGroup(grl)
  overlaps <- as.data.table(findOverlaps(grl_flanks, grl_flanks))
  overlaps <- overlaps[subjectHits > queryHits, ]
  layers <- overlaps_to_layers(overlaps)
  all_layers <- rep(1, length(grl))
  all_layers[as.numeric(names(layers))] <- layers
  all_layers <- rep(all_layers, numExonsPerGroup(grl))
  return(all_layers)
}

#'
#' @import ORFik
#' @importFrom GenomicRanges ranges resize
#' @importFrom IRanges subsetByOverlaps IRangesList
#' @keywords internal
createGeneModelPanel <- function(display_range, annotation, frame=1, custom_regions) {

  if (!is.null(custom_regions)) {
    same_names <- names(custom_regions) %in% names(annotation)
    names(custom_regions)[same_names] <- paste(names(custom_regions)[same_names], "_1", sep="")
    annotation <- c(annotation, custom_regions)
  }
  overlaps <- subsetByOverlaps(annotation, display_range)



  plot_width <- widthPerGroup(display_range)
  onames <- rep(names(overlaps), numExonsPerGroup(overlaps, FALSE))
  overlaps <- unlistGrl(overlaps)
  names(overlaps) <- onames
  overlaps$rel_frame <- getRelativeFrames(overlaps)
  overlaps <- subsetByOverlaps(overlaps, display_range)

  intersections <- trimOverlaps(overlaps,display_range)
  intersections <- groupGRangesBy(intersections)

  layers <- geneTrackLayer(intersections)



  locations <- pmapToTranscriptF(intersections, display_range)
  locations <- unlistGrl(locations)
  locations <- ranges(locations)
  blocks <- c(start(locations) , end(locations))
  names(blocks) <- rep(names(locations),2)
  blocks <- sort(blocks)
  lines_locations <- blocks[!(blocks %in% c(1, plot_width))]
  cols <- colour_bars(overlaps, display_range)
  names(cols) <- names(overlaps)
  cols <- selectCols(cols,locations)
  rect_locations <- locations

  locations <- resize(locations, width = 1, fix = "center")

  labels_locations <- start(locations)

  too_close <- labels_locations < 0.02 * plot_width

  too_far <- labels_locations > 0.98 * plot_width
  labels_locations[too_close] <- 1
  labels_locations[too_far] <- plot_width
  hjusts <- rep("center", length(labels_locations))
  hjusts[too_close] <- "left"
  hjusts[too_far] <- "right"
  if (as.logical(strand(display_range[[1]][1]) == "+")) {
    gene_names <- names(locations)
  } else gene_names <- names(locations)
  suppressWarnings({
    result_plot <- ggplot() +
      geom_rect(mapping=aes(ymin=0 - layers,ymax = 1 - layers,xmin=start(rect_locations), xmax = end(rect_locations), frame = frame),fill = cols, alpha = 0.5) +
      ylab("") +
      xlab("") +
      theme(axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank()) +
      theme(plot.margin = unit(c(0,0,0,0), "pt"))+
      scale_x_continuous(expand = c(0,0)) +
      scale_y_continuous(expand = c(0,0)) +
      geom_text(mapping = aes(y = 0.5 - layers, frame = frame, x = labels_locations, label = gene_names), color = "black", hjust = hjusts)
  })
  if (length(locations) == 0) {
    locations <- pmapToTranscriptF(display_range %>% IRangesList(),display_range)
    locations <- unlist(locations)
    result_plot <- ggplot() +
      geom_rect(mapping=aes(ymin=0,ymax=1,xmin=start(locations), xmax = end(locations)),fill = "white", alpha = 0.5) +
      ylab("") +
      xlab("") +
      theme(axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank()) +
      theme(plot.margin = unit(c(0,0,0,0), "pt"))+
      scale_x_continuous(expand = c(0,0)) +
      scale_y_continuous(expand = c(0,0))
  }
  custom_region_names <- which(names(lines_locations) %in% names(custom_regions))
  names(lines_locations) <- rep("black", length(lines_locations))
  names(lines_locations)[custom_region_names] <- "orange4"
  return(list(result_plot, lines_locations))
}

#' @importFrom Biostrings nchar
#' @keywords internal
nt_bar <- function(seq) {
  nc <- Biostrings::nchar(seq)
  position <- 1:nc
  chars <- Biostrings::strsplit(as.character(seq),"")[[1]]
  colors <- c("#619CFF","#F8766D","#00BA38")[position %% 3 + 1]
  nt_df <- data.frame(nucleotide = chars,
                      colors = colors,
                      position = position,
                      y = rep(0,nc))
  p = ggplot(nt_df, aes(x = position, y=y)) +
    geom_text(label = chars, color = colors, position = "identity") +
    # geom_point(color = colors,shape = chars) +
    theme(axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank()) +
    theme(plot.margin = unit(c(0,0,0,0), "pt"))+
    scale_x_continuous(expand = c(0,0))
  p
}


# Keep for future bars for sequence display

#' #' @importFrom Biostrings nchar
#' #' @keywords internal
#' nt_bar <- function(seq) {
#'   nc <- Biostrings::nchar(seq)
#'   position <- 1:nc
#'   chars <- Biostrings::strsplit(as.character(seq),"")[[1]]
#'   col_idx <- match(chars, c("A","T","C","G"))
#'   nt_cols <- c("yellow", "green", "red", "blue")[col_idx]
#'   colors <- c("#619CFF","#F8766D","#00BA38")[position %% 3 + 1]
#'   nt_df <- data.frame(nucleotide = chars,
#'                       colors = colors,
#'                       nt_cols = nt_cols,
#'                       start = c(1, position[-1] - 0.5),
#'                       end = c(position[-nc] + 0.5, nc),
#'                       y1 = 0,
#'                       y2 = 1)
#'   p = ggplot(nt_df) +
#'     geom_rect(aes(ymin = 0, ymax = 1, xmin = start, xmax = end), fill = nt_cols) +
#'     # geom_text(aes(y = 0.5, x = position),label = chars, color = colors, position = "identity")  +
#'     theme(axis.title = element_blank(),
#'           axis.ticks = element_blank(),
#'           axis.text = element_blank()) +
#'     theme(plot.margin = unit(c(0,0,0,0), "pt"))+
#'     scale_x_continuous(expand = c(0,0))
#'   p
#' }
