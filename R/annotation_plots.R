#'
#' @import ggplot2
#' @keywords internal
createSeqPanel <- function(sequence, start_codons = "ATG", stop_codons = c("TAA", "TAG", "TGA"), frame = 1) {
  starts <- matchMultiplePatterns(start_codons, sequence)
  stops <- matchMultiplePatterns(stop_codons, sequence)
  starts_0 <- starts[(starts - 1) %% 3 == 0 ]
  starts_1 <- starts[(starts - 1) %% 3 == 1 ]
  starts_2 <- starts[(starts - 1) %% 3 == 2 ]
  stops_0 <- stops[(stops - 1) %% 3 == 0 ]
  stops_1 <- stops[(stops - 1) %% 3 == 1 ]
  stops_2 <- stops[(stops - 1) %% 3 == 2 ]
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
  })
  return(fig)

}




#'
#' @import ORFik
#' @keywords internal
createGeneModelPanel <- function(target_range, annotation, frame=1) {
  overlaps <- subsetByOverlaps(annotation, target_range)
  plot_width <- widthPerGroup(target_range)
  overlaps <- unlistGrl(overlaps)
  overlaps$rel_frame <- getRelativeFrames(overlaps)
  overlaps <- subsetByOverlaps(overlaps, target_range)

  intersections <- trimOverlaps(overlaps,target_range)
  intersections <- groupGRangesBy(intersections)

  locations <- pmapToTranscriptF(intersections, target_range)
  locations <- unlistGrl(locations)
  locations <- ranges(locations)
  blocks <- sort(c(start(locations) , end(locations)))
  lines_locations <- blocks[!(blocks %in% c(1, plot_width))]
  cols <- colour_bars(overlaps, target_range)
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
  if (as.logical(strand(target_range[[1]][1]) == "+")) {
    gene_names <- names(locations)
  } else gene_names <- names(locations)
  suppressWarnings({
  result_plot <- ggplot() +
    geom_rect(mapping=aes(ymin=0,ymax=1,xmin=start(rect_locations), xmax = end(rect_locations), frame = frame),fill = cols, alpha = 0.5) +
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
    geom_text(mapping = aes(y=rep(0.5,length(labels_locations), frame = frame) , x = labels_locations, label = gene_names), color = "black", hjust = hjusts)
  })
  if (length(locations) == 0) {
    locations <- pmapToTranscriptF(target_range %>% IRangesList(),target_range)
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
  return(list(result_plot, lines_locations))
}
