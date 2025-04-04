#' Create sequence panel for RiboCrypt
#'
#' @param sequence the DNAStringSet
#' @param start_codons character vector, default "ATG"
#' @param stop_codons character vector, default c("TAA", "TAG", "TGA")
#' @param frame frame not used
#' @param custom_motif character vector, default NULL.
#' @return a ggplot object
#' @keywords internal
createSeqPanelPattern <- function(sequence, start_codons = "ATG", stop_codons = c("TAA", "TAG", "TGA"), frame = 1,
                           custom_motif = NULL) {
  if (is.null(custom_motif)) {
  } else if (custom_motif == "") {
    custom_motif <- NULL
  } else {
    custom_motif <- toupper(custom_motif)
    custom_motif <- gsub("U", "T", custom_motif)
    custom_motif <- strsplit(custom_motif, ",")[[1]] %>% strsplit(" ") %>% unlist
  }
  hits <- lapply(list(start_codons, stop_codons, custom_motif), function(x) matchMultiplePatterns(x, sequence))
  names(hits) <- c("white", "black", "#C97D00")
  hits <- lapply(hits, as.data.table)
  hits <- rbindlist(hits, idcol = "col")
  names(hits) <- c("col", "pos")
  pos <- NULL # avoid dt warning
  hits[, frames := (pos - 1) %% 3]

  return(hits)

}

plotSeqPanel <- function(hits, sequence, frame = 1) {
  pos <- NULL # avoid dt warning
  # New red: #FFA8A3 ?
  fig <- ggplot() +
    geom_rect(aes(ymin = c(1,0,-1), ymax = c(2,1,0), xmin = rep(1,3), xmax = rep(length(sequence),3)),
              fill = c("#F8766D","#00BA38","#619CFF")) +
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
#'
#' @importFrom IRanges findOverlaps
#' @param grl a GRangesList
#' @return numeric, the track row index
#' @keywords internal
geneTrackLayer <- function(grl) {
  grl_flanks <- flankPerGroup(grl)
  overlaps <- as.data.table(findOverlaps(grl_flanks, grl_flanks))
  overlaps <- overlaps[subjectHits > queryHits, ]

  if (nrow(overlaps) == 0) {
    all_layers <- rep(1, length(grl))
    all_layers <- rep(all_layers, numExonsPerGroup(grl))
  } else {
    layers <- overlaps_to_layers(overlaps)
    all_layers <- rep(1, length(grl))
    all_layers[as.numeric(names(layers))] <- layers
    all_layers <- rep(all_layers, numExonsPerGroup(grl))
  }
  return(all_layers)
}

#' Create gene annotation layer for browser
#' @import ORFik
#' @importFrom GenomicRanges ranges resize psetdiff
#' @importFrom IRanges subsetByOverlaps IRangesList
#' @keywords internal
#' @noRd
createGeneModelPanel <- function(display_range, annotation, tx_annotation = NULL,
                                 frame=1, custom_regions, viewMode) {
  # TODO: Explain sections of this function, or split in sub functions
  # It is too complicated right now.
  use_custom_region <- !is.null(custom_regions) & length(custom_regions) > 0
  if (use_custom_region) {
    same_names <- names(custom_regions) %in% names(annotation)
    names(custom_regions)[same_names] <- paste(names(custom_regions)[same_names], "_1", sep="")
    overlaps_custom <- subsetByOverlaps(custom_regions, display_range,
                                        type = ifelse(viewMode == "tx", "within", "any"))
  }
  overlaps <- subsetByOverlaps(annotation, display_range,
                               type = ifelse(viewMode == "tx", "within", "any"))
  overlaps_tx <- NULL
  if (viewMode != "tx") {
    overlaps_tx <- subsetByOverlaps(tx_annotation, display_range,
                                    type = ifelse(viewMode == "tx", "within", "any"))
    if (length(overlaps) > 0) {
      if (!all(names(overlaps_tx) %in% names(overlaps))) {
        overlaps_tx <- overlaps_tx[names(overlaps_tx) %in% names(overlaps)]
      }
      overlaps_tx <- groupGRangesBy(unlistGrl(GenomicRanges::psetdiff(unlistGrl(overlaps_tx),
                                                                      overlaps[names(unlistGrl(overlaps_tx))])))
    }
  }

  if (use_custom_region) overlaps <- c(overlaps, overlaps_custom)
  if (length(overlaps) > 0) overlaps@unlistData$type <- "cds"
  if (length(overlaps_tx) > 0) overlaps_tx@unlistData$type <- "utr"

  return(gene_box_fix_overlaps(display_range, c(overlaps, overlaps_tx), custom_regions))
}

gene_box_fix_overlaps <- function(display_range, overlaps, custom_regions) {
  if (length(overlaps) == 0) {
    result_dt <- data.table()
    lines_locations <- NULL
    return(list(result_dt, lines_locations))
  }

  names_grouping <- rep(names(overlaps), lengths(overlaps))
  overlaps <- unlistGrl(overlaps)
  names(overlaps) <- names_grouping
  overlaps$rel_frame_exon <- getRelativeFrames(overlaps)
  overlaps <- subsetByOverlaps(overlaps, display_range)


  intersections <- trimOverlaps(overlaps, display_range)
  intersections <- groupGRangesBy(intersections, paste0(names(intersections), "___", intersections$type))
  names(intersections) <- sub("___.*", "", names(intersections))

  locations <- pmapToTranscriptF(intersections, display_range)
  type_per_grl <- unlist(lapply(intersections, function(x) unique(x$type)))
  locations@unlistData$type <- rep(type_per_grl, times = lengths(locations))
  locations <- unlistGrl(locations)
  locations$rel_frame_orf <- getRelativeFrames(locations)
  locations <- unlistGrl(groupGRangesBy(locations))

  cols <- colour_bars(locations, overlaps, display_range)
  layers <- geneTrackLayer(groupGRangesBy(locations))
  plot_width <- widthPerGroup(display_range)
  geneBox <- geneBoxFromRanges(locations, plot_width, layers, cols,
                               custom_regions)
  return(geneBox)
}

#' Get gene annotation box coordinates
#'
#' Get data.table of coordinates and coloring, and a lines vector for
#' main plot vertical lines.
#' @noRd
geneBoxFromRanges <- function(locations, plot_width,
                              layers = rep(1, length(locations)),
                              cols = rc_rgb()[start(locations) %% 3 + 1],
                              custom_regions = NULL) {
  type <- locations$type
  locations <- move_utrs_touching_cds(locations, type, plot_width)


  locations <- ranges(locations)
  blocks <- c(start(locations) , end(locations))
  names(blocks) <- rep(names(locations), 2)
  blocks <- sort(blocks)
  lines_locations <- blocks[!(blocks %in% c(1, plot_width))]


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
                          gene_names = gene_names,
                          type = type)
  result_dt[,no_ex := .N, by = gene_names]
  result_dt <- result_dt[order(type, decreasing = TRUE),]
  dt <- result_dt
  draw_introns <- nrow(dt[no_ex > 1]) > 0
  if (draw_introns) {
    seg_dt <- dt[no_ex > 1]
    seg_dt <- seg_dt[ , .(layers = layers[1], rect_starts = rect_ends[1:(.N - 1)], rect_ends = rect_starts[2:.N],
                          cols = "grey45", labels_locations = 0, hjusts = "center",
                          type = "intron", no_ex = "1"), by = gene_names]
    seg_dt <- seg_dt[, colnames(dt), with = FALSE]
    with_introns_dt <- rbindlist(list(dt, seg_dt))
    result_dt <- with_introns_dt
  }
  res_list <- list(result_dt, lines_locations)
  return(res_list)
}

move_utrs_touching_cds <- function(locations, type, plot_width) {
  utr_end_touch_cds <- type == "utr" & (end(locations)+1) %in% start(locations)[type == "cds"]
  new_utr_ends <- end(locations)[utr_end_touch_cds] + 1
  end(locations)[utr_end_touch_cds] <- pmin(plot_width, new_utr_ends)

  utr_start_touch_cds <- type == "utr" & (start(locations)-1) %in% end(locations)[type == "cds"]
  new_utr_starts <- start(locations)[utr_start_touch_cds] - 1
  start(locations)[utr_start_touch_cds] <- pmax(1, new_utr_starts)
  return(locations)
}

rc_rgb <- function() {
  return(c("#fbbab5","#7fdc9b","#afcdff"))
}

geneModelPanelPlot <- function(dt, frame = 1) {
  # If no annotation given, dt is empty, return blank
  if (nrow(dt) == 0) return(ggplot())
  # Else render exon boxes
  seg_dt <- dt[type == "intron"]
  dt <- dt[type != "intron"]

  result_plot <- ggplot(frame = frame) + ylab("") + xlab("") +
    theme(axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.margin = unit(c(0,0,0,0), "pt")) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    theme(panel.background = element_rect(fill= "white"))

  draw_introns <- nrow(seg_dt) > 0
  if (draw_introns) result_plot <- result_plot +
    geom_segment(data = seg_dt,
                 mapping = aes(x = rect_starts, xend = rect_ends, y = 0.5 - layers, yend = 0.5 - layers, text = gene_names),
                 color = "grey45", alpha = 0.6)

  suppressWarnings({
    result_plot <-  result_plot +
     geom_rect(data = dt, mapping=aes(ymin=0 - layers + ifelse(type == "cds", 0, 0.33), ymax = 1 - layers - ifelse(type == "cds", 0, 0.33),
                                      xmin=rect_starts,xmax = rect_ends, text = gene_names),
               fill = dt$cols, color = "grey45") +
      geom_text(data = dt[,.(layers = layers[1], labels_locations = mean(labels_locations)),gene_names],
                mapping = aes(y = 0.50 - layers, x = labels_locations,
                              label = gene_names), color = "black", hjust = "center", vjust = "center")
    })

  return(result_plot)
}

nt_area_template <- function() {
  ggplot() +
    theme(axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank()) +
    theme(plot.margin = unit(c(0,0,0,0), "pt"))+
    scale_x_continuous(expand = c(0,0))
}
