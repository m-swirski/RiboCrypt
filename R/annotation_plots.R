#' Create sequence panel for RiboCrypt
#'
#' @param start_codons character vector, default "ATG"
#' @param stop_codons character vector, default c("TAA", "TAG", "TGA")
#' @param custom_motif character vector, default NULL.
#' @import ggplot2
#' @return a ggplot object
#' @keywords internal
createSeqPanelPattern <- function(sequence, start_codons = "ATG", stop_codons = c("TAA", "TAG", "TGA"), frame = 1,
                           custom_motif = NULL) {

  hits <- lapply(list(start_codons, stop_codons, custom_motif), function(x) matchMultiplePatterns(x, sequence))
  names(hits) <- c("white", "black", "purple")
  hits <- lapply(hits, as.data.table)
  hits <- rbindlist(hits, idcol = "col")
  names(hits) <- c("col", "pos")
  hits[, frames := pos %% 3]

  return(hits)

}

plotSeqPanel <- function(hits, sequence, frame = 1) {
  fig <- ggplot() +
    geom_rect(aes(ymin = c(1,0,-1), ymax = c(2,1,0), xmin = rep(1,3), xmax = rep(length(sequence),3), frame = frame),
              fill = c("#F8766D","#00BA38","#619CFF" )) +
    geom_segment(data=hits, mapping = aes(y = 2 - (frames + 1), yend =  2 - frames, x = pos, xend = pos), col=hits$col) +
    ylab("frame") +
    xlab("position [nt]") +
    theme_bw() +
    theme(plot.margin = unit(c(0,0,0,0), "pt"), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
    scale_y_continuous(breaks = c(-0.5,0.5, 1.5), labels = c("2","1", "0"), expand = c(0,0)) + scale_x_continuous(expand = c(0,0))

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

#' How many rows does the gene track need
#' @keywords internal
geneTrackLayer <- function(grl) {

  grl_flanks <- flankPerGroup(grl)
  overlaps <- as.data.table(findOverlaps(grl_flanks, grl_flanks))
  overlaps <- overlaps[subjectHits > queryHits, ]

  if (nrow(overlaps) == 0) {
    all_layers <- rep(1, length(grl))
    all_layers <- rep(all_layers, numExonsPerGroup(grl))
    return(all_layers)
  } else {
  layers <- overlaps_to_layers(overlaps)
  all_layers <- rep(1, length(grl))
  all_layers[as.numeric(names(layers))] <- layers
  all_layers <- rep(all_layers, numExonsPerGroup(grl))
  return(all_layers)
  }
}

#'
#' @import ORFik
#' @importFrom GenomicRanges ranges resize
#' @importFrom IRanges subsetByOverlaps IRangesList
#' @keywords internal


createGeneModelPanel <- function(display_range, annotation, frame=1, custom_regions, viewMode) {
# TODO: Explain sections of this function, or split in sub functions
# It is too complicated right now.
if (!is.null(custom_regions)) {
  same_names <- names(custom_regions) %in% names(annotation)
  names(custom_regions)[same_names] <- paste(names(custom_regions)[same_names], "_1", sep="")
  annotation <- c(annotation, custom_regions)
}
overlaps <- subsetByOverlaps(annotation, display_range,
                             type = ifelse(viewMode == "tx", "within", "any"))
# Default theme


if (length(overlaps) > 0) {

  plot_width <- widthPerGroup(display_range)
  onames <- rep(names(overlaps), numExonsPerGroup(overlaps, FALSE))
  overlaps <- unlistGrl(overlaps)
  names(overlaps) <- onames
  overlaps$rel_frame <- getRelativeFrames(overlaps)
  overlaps <- subsetByOverlaps(overlaps, display_range)

  intersections <- trimOverlaps(overlaps,display_range)
  intersections <- groupGRangesBy(intersections)

  locations <- pmapToTranscriptF(intersections, display_range)
  layers <- geneTrackLayer(locations)

  locations <- unlistGrl(locations)
  locations <- ranges(locations)
  blocks <- c(start(locations) , end(locations))
  names(blocks) <- rep(names(locations),2)
  blocks <- sort(blocks)
  lines_locations <- blocks[!(blocks %in% c(1, plot_width))]
  cols <- colour_bars(overlaps, display_range)
  names(cols) <- names(overlaps)
  if (length(cols) != length(locations)) cols <- selectCols(cols,locations)
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
  # TODO: remove when verified this is not needed
  # if (as.logical(strand(display_range[[1]][1]) == "+")) {
  #   gene_names <- names(locations)
  # } else gene_names <- names(locations)
  gene_names <- names(locations)
  custom_region_names <- which(names(lines_locations) %in% names(custom_regions))
  names(lines_locations) <- rep("black", length(lines_locations))
  names(lines_locations)[custom_region_names] <- "orange4"
  result_dt <- data.table(layers = layers,
                          rect_starts = start(rect_locations),
                          rect_ends = end(rect_locations),
                          cols = cols,
                          labels_locations = labels_locations,
                          hjusts = hjusts,
                          gene_names = gene_names)

} else {
  result_dt <- data.table()
  lines_locations <- NULL
}

return(list(result_dt, lines_locations))
}

geneModelPanelPlot <- function(dt, frame = 1) {

  base_gg <- ggplot(frame = frame) + ylab("") + xlab("") +
    theme(axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.margin = unit(c(0,0,0,0), "pt")) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    theme(panel.background = element_rect(fill= "white"))

  if (nrow(dt) > 0) {
    suppressWarnings({
      result_plot <- base_gg +
        geom_rect(data = dt, mapping=aes(ymin=0 - layers, ymax = 1 - layers, xmin=rect_starts,
                                         xmax = rect_ends),
                  fill = dt$cols, alpha = 0.5) +
        geom_text(data = dt, mapping = aes(y = 0.5 - layers, x = labels_locations,
                                           label = gene_names), color = "black", hjust = dt$hjusts)})
  } else result_plot <- base_gg
  return(result_plot)
}
#'
#' createGeneModelPanel <- function(display_range, annotation, frame=1, custom_regions, viewMode) {
#'   # TODO: Explain sections of this function, or split in sub functions
#'   # It is too complicated right now.
#'   if (!is.null(custom_regions)) {
#'     same_names <- names(custom_regions) %in% names(annotation)
#'     names(custom_regions)[same_names] <- paste(names(custom_regions)[same_names], "_1", sep="")
#'     annotation <- c(annotation, custom_regions)
#'   }
#'   overlaps <- subsetByOverlaps(annotation, display_range,
#'                                type = ifelse(viewMode == "tx", "within", "any"))
#'   # Default theme
#'   base_gg <- ggplot(frame = frame) + ylab("") + xlab("") +
#'     theme(axis.title.x = element_blank(),
#'           axis.ticks.x = element_blank(),
#'           axis.text.x = element_blank(),
#'           axis.text.y = element_blank(),
#'           axis.ticks.y = element_blank(),
#'           plot.margin = unit(c(0,0,0,0), "pt")) +
#'     scale_x_continuous(expand = c(0,0)) +
#'     scale_y_continuous(expand = c(0,0))
#'
#'   if (length(overlaps) > 0) {
#'
#'   plot_width <- widthPerGroup(display_range)
#'   onames <- rep(names(overlaps), numExonsPerGroup(overlaps, FALSE))
#'   overlaps <- unlistGrl(overlaps)
#'   names(overlaps) <- onames
#'   overlaps$rel_frame <- getRelativeFrames(overlaps)
#'   overlaps <- subsetByOverlaps(overlaps, display_range)
#'
#'   intersections <- trimOverlaps(overlaps,display_range)
#'   intersections <- groupGRangesBy(intersections)
#'
#'   locations <- pmapToTranscriptF(intersections, display_range)
#'   layers <- geneTrackLayer(locations)
#'
#'   locations <- unlistGrl(locations)
#'   locations <- ranges(locations)
#'   blocks <- c(start(locations) , end(locations))
#'   names(blocks) <- rep(names(locations),2)
#'   blocks <- sort(blocks)
#'   lines_locations <- blocks[!(blocks %in% c(1, plot_width))]
#'   cols <- colour_bars(overlaps, display_range)
#'   names(cols) <- names(overlaps)
#'   cols <- selectCols(cols,locations)
#'   rect_locations <- locations
#'
#'   locations <- resize(locations, width = 1, fix = "center")
#'
#'   labels_locations <- start(locations)
#'
#'   too_close <- labels_locations < 0.02 * plot_width
#'
#'   too_far <- labels_locations > 0.98 * plot_width
#'   labels_locations[too_close] <- 1
#'   labels_locations[too_far] <- plot_width
#'   hjusts <- rep("center", length(labels_locations))
#'   hjusts[too_close] <- "left"
#'   hjusts[too_far] <- "right"
#'   # TODO: remove when verified this is not needed
#'   # if (as.logical(strand(display_range[[1]][1]) == "+")) {
#'   #   gene_names <- names(locations)
#'   # } else gene_names <- names(locations)
#'   gene_names <- names(locations)
#'   custom_region_names <- which(names(lines_locations) %in% names(custom_regions))
#'   names(lines_locations) <- rep("black", length(lines_locations))
#'   names(lines_locations)[custom_region_names] <- "orange4"
#'   suppressWarnings({
#'     result_plot <- base_gg +
#'       geom_rect(mapping=aes(ymin=0 - layers, ymax = 1 - layers, xmin=start(rect_locations),
#'                             xmax = end(rect_locations)),
#'                 fill = cols, alpha = 0.5) +
#'       geom_text(mapping = aes(y = 0.5 - layers, x = labels_locations,
#'                               label = gene_names), color = "black", hjust = hjusts)
#'   })
#'   } else {
#'     locations <- pmapToTranscriptF(display_range %>% IRangesList(),display_range)
#'     locations <- unlist(locations)
#'     result_plot <- base_gg +
#'       geom_rect(mapping=aes(ymin=0,ymax=1,xmin=start(locations), xmax = end(locations)),
#'                 fill = "white", alpha = 0.5)
#'     lines_locations <- NULL
#'   }
#'
#'   return(list(result_plot, lines_locations))
#' }
#'
#' #' @importFrom Biostrings nchar
#' #' @keywords internal
#' nt_bar <- function(seq) {
#'   nc <- Biostrings::nchar(seq)
#'   position <- 1:nc
#'   chars <- Biostrings::strsplit(as.character(seq),"")[[1]]
#'   colors <- c("#619CFF","#F8766D","#00BA38")[position %% 3 + 1]
#'   nt_df <- data.frame(nucleotide = chars,
#'                       colors = colors,
#'                       position = position,
#'                       y = rep(0,nc))
#'   p = ggplot(nt_df, aes(x = position, y=y)) +
#'     geom_text(label = chars, color = colors, position = "identity") +
#'     # geom_point(color = colors,shape = chars) +
#'     theme(axis.title = element_blank(),
#'           axis.ticks = element_blank(),
#'           axis.text = element_blank()) +
#'     theme(plot.margin = unit(c(0,0,0,0), "pt"))+
#'     scale_x_continuous(expand = c(0,0))
#'   p
#' }
#'

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
