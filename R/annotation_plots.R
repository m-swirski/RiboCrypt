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

plotAASeqPanel <- function(hits, sequence, frame_colors = "R", gg_theme = theme_bw()) {
  stopifnot(is(gg_theme, "theme"))
  frame_colors <- frame_color_themes(frame_colors, FALSE)
  pos <- NULL # avoid dt warning
  # New red: #FFA8A3 ?
  hits[, Type := fifelse(col == "white", "Start codon", fifelse(col == "black", "Stop codon", "User Motif"))]
  seq_panel_template <- attr(gg_theme, "seq_panel_template")

  if (is.null(seq_panel_template)) seq_panel_template <- seqPanelPlotTemplate()
  # theme_bw propogates to all subplots
  fig <- seq_panel_template +
    geom_rect(aes(ymin = c(1,0,-1), ymax = c(2,1,0), xmin = rep(1,3), xmax = rep(length(sequence),3)),
              fill = frame_colors) +
    suppressWarnings(
      geom_segment(data=hits,
                   mapping = aes(y = 2 - (frames + 1), yend =  2 - frames, x = pos, xend = pos,
                                 text = paste0("Pos: ", pos, "\n", "Type: ", Type)), col=hits$col)
      )

  return(fig)
}

gg_theme_template <- function() {
  gg_theme <- theme_bw()
  attr(gg_theme, "seq_panel_template") <- seqPanelPlotTemplate(gg_theme)
  attr(gg_theme, "gg_template") <- geneModelPanelPlotTemplate()
  attr(gg_theme, "seq_panel_nt_template_ggplot") <- nt_area_template()
  attr(gg_theme, "seq_panel_nt_template_plotly") <-
    automateTicksLetters(attr(gg_theme, "seq_panel_nt_template_ggplot"))
  return(gg_theme)
}

seqPanelPlotTemplate <- function(gg_theme) {
  suppressWarnings(ggplot(frame = "static")) +
    ylab("frame") +
    xlab("position [nt]") +
    gg_theme +
    theme(plot.margin = unit(c(0,0,0,0), "pt"), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
    scale_y_continuous(breaks = c(-0.5,0.5, 1.5), labels = c("2","1", "0"), expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0))
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
                                 custom_regions, viewMode, collapse_intron_flank = 100,
                                 frame_colors = "R") {
  # TODO: Explain sections of this function, or split in sub functions
  # It is too complicated right now.
  use_custom_region <- !is.null(custom_regions) & length(custom_regions) > 0
  overlap_type <- ifelse(viewMode == "tx", "within", "any")
  if (use_custom_region) {
    same_names <- names(custom_regions) %in% names(annotation)
    names(custom_regions)[same_names] <- paste(names(custom_regions)[same_names], "_1", sep="")
    overlaps_custom <- subsetByOverlaps(custom_regions, display_range,
                                        type = overlap_type)
  }

  if (all(length(annotation) == 1 & length(display_range) == 1) & identical(names(annotation), names(display_range))) {
    annotation_shown <- annotation
  } else annotation_shown <- subsetByOverlaps(annotation, display_range,
                                      type = overlap_type)

  overlaps_tx <- NULL

  if (viewMode != "tx") {
    if (all(length(tx_annotation) == 1 & length(display_range) == 1) & identical(names(tx_annotation), names(display_range))) {
      overlaps_tx <- tx_annotation
    } else overlaps_tx <- subsetByOverlaps(tx_annotation, display_range,
                                        type = overlap_type)
    if (length(annotation_shown) > 0) {
      if (!all(names(overlaps_tx) %in% names(annotation_shown))) {
        overlaps_tx <- overlaps_tx[names(overlaps_tx) %in% names(annotation_shown)]
      }
      unl_grl <- unlistGrl(overlaps_tx)
      overlaps_tx <- groupGRangesBy(
        unlistGrl(GenomicRanges::psetdiff(unl_grl, annotation_shown[names(unl_grl)])))
    }
  }

  if (use_custom_region) annotation_shown <- c(annotation_shown, overlaps_custom)
  if (length(annotation_shown) > 0) annotation_shown@unlistData$type <- "cds"
  if (length(overlaps_tx) > 0) overlaps_tx@unlistData$type <- "utr"

  return(gene_box_fix_overlaps(display_range, c(annotation_shown, overlaps_tx), custom_regions,
                               collapse_intron_flank, frame_colors))
}

gene_box_fix_overlaps <- function(display_range, overlaps, custom_regions,
                                  collapse_intron_flank, frame_colors) {
  if (length(overlaps) == 0) {
    result_dt <- data.table()
    lines_locations <- NULL
    return(list(result_dt, lines_locations))
  }

  names_grouping <- rep(names(overlaps), lengths(overlaps))
  overlaps <- unlistGrl(overlaps)
  names(overlaps) <- names_grouping
  overlaps$rel_frame_exon <- getRelativeFrames(overlaps)
  # overlaps <- subsetByOverlaps(overlaps, display_range)

  intersections <- trimOverlaps(overlaps, display_range)
  intersections <- groupGRangesBy(intersections, paste0(names(intersections), "___", intersections$type))
  names(intersections) <- sub("___.*", "", names(intersections))

  locations <- pmapToTranscriptF(intersections, display_range)
  type_per_grl <- unlist(lapply(intersections, function(x) unique(x$type)))
  locations@unlistData$type <- rep(type_per_grl, times = lengths(locations))
  locations <- unlistGrl(locations)
  locations$rel_frame_orf <- getRelativeFrames(locations)
  locations <- unlistGrl(groupGRangesBy(locations))

  cols <- colour_bars(locations, overlaps, display_range, frame_colors)
  layers <- geneTrackLayer(groupGRangesBy(locations))
  plot_width <- widthPerGroup(display_range)
  geneBox <- geneBoxFromRanges(locations, plot_width, layers, cols,
                               custom_regions, collapse_intron_flank)
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
                              custom_regions = NULL, collapse_intron_flank = 100) {
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
  gene_names <- names(locations)

  cap <- 0.00133*nchar(gene_names)
  too_close <- labels_locations < cap * plot_width
  too_far <- labels_locations > (1-cap) * plot_width

  labels_locations[too_close] <- 1
  labels_locations[too_far] <- plot_width
  hjusts <- rep("center", length(labels_locations))
  hjusts[too_close] <- "left"
  hjusts[too_far] <- "right"

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
    ranges <- reduce(GRanges(seg_dt$gene_names, ranges = IRanges(seg_dt$rect_starts, seg_dt$rect_ends)))
    did_collapse_introns <- (nrow(seg_dt) > 1 & length(ranges) != 1) & collapse_intron_flank > 0
    if (did_collapse_introns ) {
      seg_dt <- data.table(gene_names = as.character(seqnames(ranges)), rect_starts = start(ranges), rect_ends = end(ranges))
      unique_layers <- dt[, .(layers = unique(layers)[1]), by = gene_names]
      seg_dt <- merge.data.table(seg_dt, unique_layers, by = "gene_names", sort = FALSE)
    }

    seg_dt <- seg_dt[ , .(layers = layers[1], rect_starts = rect_ends[1:(.N - 1)], rect_ends = rect_starts[2:.N],
                          cols = "grey45", labels_locations = 0, hjusts = "center",
                          type = "intron", no_ex = "1"), by = gene_names]
    if (did_collapse_introns) {
      dt_unique <- dt[!duplicated(dt[, .(rect_starts, rect_ends)]), ][, .(rect_starts, rect_ends)]
      dt_unique_all_pos <- unique(unlist(lapply(as.data.table(t(as.matrix(dt_unique))), function(x) {seq.int(x[1], x[2])}), use.names = FALSE))
      seg_dt_all_pos <- lapply(as.data.table(t(as.matrix(seg_dt[, .(rect_starts, rect_ends)]))), function(x) {seq.int(x[1]+1, x[2]-1)})
      overlaps_other_exon <- sapply(seq_along(seg_dt_all_pos), function(i) {any(seg_dt_all_pos[[i]] %in% dt_unique_all_pos)})
      seg_dt[(rect_ends - rect_starts > collapse_intron_flank*2) & !overlaps_other_exon, type := "intron_collapsed"]
    }
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

frame_color_themes <- function(theme, with_alpha = FALSE) {
  if (theme == "R") {
    colors <- rc_rgb(with_alpha)
  } else if (theme == "Color_blind") {
    colors <- if (with_alpha) {
      c("#8B90BA", "#FBCF79", "#C6E0EE")
    } else c("#5056A0", "#F9AA2A", "#AED5EB")
  } else stop("theme must be either R or Color blind")
}

rc_rgb <- function(with_alpha = TRUE) {
  rgb_colors <- if (with_alpha) {
    c("#fbbab5","#7fdc9b","#afcdff")
  } else c("#F8766D","#00BA38","#619CFF")
  return(rgb_colors)
}

geneModelPanelPlot <- function(dt, gg_template = geneModelPanelPlotTemplate()) {
  # If no annotation given, dt is empty, return blank
  if (nrow(dt) == 0) return(ggplot())
  # Else render exon boxes
  seg_dt <- dt[type %in% c("intron", "intron_collapsed")]
  dt <- dt[!(type %in% c("intron", "intron_collapsed"))]

  result_plot <- gg_template

  trans <- c("cds", "translon")
  draw_introns <- nrow(seg_dt) > 0
  if (draw_introns)  {
    intron_flank_coords <- seg_dt[type %in% c("intron_collapsed")]
    intron_flank_coords[, tooltip := "Collapsed Intron"]
    suppressWarnings({
    result_plot <- result_plot +
      geom_segment(data = seg_dt,
                   mapping = aes(x = rect_starts, xend = rect_ends, y = 0.5 - layers, yend = 0.5 - layers, text = gene_names),
                   color = "grey45", alpha = 0.6) +
      geom_text(data = intron_flank_coords, aes(x = (rect_starts + rect_ends) / 2, y = 0.5 - layers, text = tooltip),
                label = "...", size = 6)
    })
  }

  suppressWarnings({
    result_plot <-  result_plot +
     geom_rect(data = dt, mapping=aes(ymin=0 - layers + ifelse(type %in% trans, 0, 0.33), ymax = 1 - layers - ifelse(type %in% trans, 0, 0.33),
                                      xmin=rect_starts,xmax = rect_ends, text = gene_names),
               fill = dt$cols, color = "grey45") +
      geom_text(data = dt[,.(layers = layers[1], labels_locations = mean(labels_locations)), gene_names],
                mapping = aes(y = 0.50 - layers, x = labels_locations,
                              label = gene_names), color = "black", hjust = "center", vjust = "center")
    })

  return(result_plot)
}

geneModelPanelPlotTemplate <- function() {
  suppressWarnings(ggplot(frame = "static")) +
    ylab("") + xlab("") +
    theme(axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.margin = unit(c(0,0,0,0), "pt")) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    theme(panel.background = element_rect(fill= "white"))
}

# Pure plotly version (no ggplot/ggplotly)
geneModelPanelPlotly <- function(dt) {
  # If no annotation given, return blank plotly
  if (is.null(dt) || nrow(dt) == 0) {
    return(plotly::plot_ly() |>
             plotly::layout(
               xaxis = list(visible = FALSE),
               yaxis = list(visible = FALSE),
               margin = list(l = 0, r = 0, b = 0, t = 0),
               paper_bgcolor = "white",
               plot_bgcolor  = "white"
             ))
  }

  stopifnot(is(dt, "data.table"))

  # Split introns vs boxes
  seg_dt <- dt[type %in% c("intron", "intron_collapsed")]
  box_dt <- dt[!(type %in% c("intron", "intron_collapsed"))]

  trans <- c("cds", "translon")

  # Compute rectangle y coords to match ggplot logic:
  # ymin = 0 - layers + ifelse(type %in% trans, 0, 0.33)
  # ymax = 1 - layers - ifelse(type %in% trans, 0, 0.33)
  box_dt[, `:=`(
    ymin = 0 - layers + ifelse(type %in% trans, 0, 0.33),
    ymax = 1 - layers - ifelse(type %in% trans, 0, 0.33)
  )]

  # Helpful for axis range
  x_min <- min(dt$rect_starts, dt$rect_ends, na.rm = TRUE)
  x_max <- max(dt$rect_starts, dt$rect_ends, na.rm = TRUE)
  y_min <- min(box_dt$ymin, if (nrow(seg_dt)) (0.5 - seg_dt$layers) else Inf, na.rm = TRUE) - 0.05
  y_max <- max(box_dt$ymax, if (nrow(seg_dt)) (0.5 - seg_dt$layers) else -Inf, na.rm = TRUE) + 0.05

  p <- plotly::plot_ly() |>
    plotly::layout(
      showlegend = FALSE,
      xaxis = list(
        visible = FALSE,
        range = c(x_min, x_max),
        fixedrange = FALSE,
        zeroline = FALSE
      ),
      yaxis = list(
        visible = FALSE,
        range = c(y_min, y_max),
        fixedrange = TRUE,
        zeroline = FALSE
      ),
      margin = list(l = 0, r = 0, b = 0, t = 0),
      paper_bgcolor = "white",
      plot_bgcolor  = "white",
      hovermode = "closest"
    )

  # --- Introns as segments (grey lines) ---
  if (nrow(seg_dt) > 0) {
    # One trace for all segments
    seg_x <- c(rbind(seg_dt$rect_starts, seg_dt$rect_ends, NA_real_))
    seg_y <- c(rbind(0.5 - seg_dt$layers, 0.5 - seg_dt$layers, NA_real_))

    # Hover text per segment: repeat text for start/end, NA for breaks
    seg_text <- c(rbind(seg_dt$gene_names, seg_dt$gene_names, NA_character_))

    p <- p |>
      plotly::add_trace(
        x = seg_x, y = seg_y,
        type = "scatter", mode = "lines",
        line = list(color = "grey45", width = 2),
        opacity = 0.6,
        text = seg_text,
        hoverinfo = "text"
      )

    # Collapsed intron label "..."
    intron_flank_coords <- seg_dt[type %in% c("intron_collapsed")]
    if (nrow(intron_flank_coords) > 0) {
      intron_flank_coords[, `:=`(
        xmid = (rect_starts + rect_ends) / 2,
        ymid = 0.5 - layers,
        tooltip = "Collapsed Intron"
      )]

      p <- p |>
        plotly::add_text(
          data = intron_flank_coords,
          x = ~xmid, y = ~ymid,
          text = "...",
          textfont = list(size = 18, color = "black"),
          hovertext = ~tooltip,
          hoverinfo = "text"
        )
    }
  }

  # --- Exon/CDS/etc boxes as filled rectangles ---
  if (nrow(box_dt) > 0) {
    # Shapes for rectangles (fast render, but shapes don't carry per-shape hover)
    rect_shapes <- lapply(seq_len(nrow(box_dt)), function(i) {
      list(
        type = "rect",
        layer = "below",
        xref = "x", yref = "y",
        x0 = box_dt$rect_starts[i], x1 = box_dt$rect_ends[i],
        y0 = box_dt$ymin[i],        y1 = box_dt$ymax[i],
        line = list(color = "grey45", width = 1),
        fillcolor = box_dt$cols[i],
        opacity = 1
      )
    })

    p <- p |>
      plotly::layout(shapes = rect_shapes)

    # Add an invisible scatter trace to provide per-rectangle hover tooltips
    box_dt[, `:=`(
      xmid = (rect_starts + rect_ends) / 2,
      ymid = (ymin + ymax) / 2
    )]

    p <- p |>
      plotly::add_markers(
        data = box_dt,
        x = ~xmid, y = ~ymid,
        marker = list(opacity = 0),
        hovertext = ~gene_names,
        hoverinfo = "text"
      )

    # Gene name labels per (layer, gene_names) at mean(labels_locations)
    lab_dt <- box_dt[, .(
      layers = layers[1],
      labels_locations = mean(labels_locations)
    ), by = gene_names]

    p <- p |>
      plotly::add_text(
        data = lab_dt,
        x = ~labels_locations,
        y = ~(0.50 - layers),
        text = ~gene_names,
        textposition = "middle center",
        textfont = list(color = "black", size = 16),
        hoverinfo = "skip"
      )
  }
  p %>% plotly::config(displayModeBar = FALSE) %>%
    plotly::layout(margin = margin_megabrowser())
}

margin_megabrowser <- function() {
  list(
    l = 30,
    r = 100,
    t = 0,
    b = 0
  )
}


nt_area_template <- function() {
  ggplot() +
    theme(axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank()) +
    theme(plot.margin = unit(c(0,0,0,0), "pt"))+
    scale_x_continuous(expand = c(0,0))
}
