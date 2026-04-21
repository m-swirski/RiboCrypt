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
    custom_motif <- gsub("U", "T", custom_motif, fixed = TRUE)
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

plotAASeqPanel <- function(hits, sequence, frame_colors = "R", theme_template = theme_bw()) {
  stopifnot(is(theme_template, "theme"))
  frame_colors <- frame_color_themes(frame_colors, FALSE)
  pos <- NULL # avoid dt warning
  # New red: #FFA8A3 ?
  hits[, Type := fifelse(col == "white", "Start codon", fifelse(col == "black", "Stop codon", "User Motif"))]
  seq_panel_template <- attr(theme_template, "seq_panel_template")

  if (is.null(seq_panel_template)) seq_panel_template <- seqPanelPlotTemplate(theme_template)
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

plotAASeqPanelPlotly <- function(hits, sequence, frame_colors = "R", template = NULL) {
  seq_length <- nchar(sequence)
  if (inherits(template, "plotly")) {
    p <- template
  } else {
    p <- aaSeqPanelPlotlyTemplate(frame_colors = frame_colors)
  }
  p <- aaSeqPanelPlotlySetFrameColors(p, frame_colors)
  p <- aaSeqPanelPlotlySetRange(p, seq_length)

  if (nrow(hits) > 0) {
    hits <- data.table::copy(hits)
    hits[, Type := fifelse(col == "white", "Start codon",
                           fifelse(col == "black", "Stop codon", "User Motif"))]
    hits[, `:=`(
      y = 1 - frames,
      yend = 2 - frames,
      hover_text = paste0("Pos: ", pos, "\nType: ", Type)
    )]

    for (hit_col in unique(hits$col)) {
      hit_dt <- hits[col == hit_col]
      p <- plotly::add_segments(
        p,
        inherit = FALSE,
        data = hit_dt,
        x = ~pos,
        xend = ~pos,
        y = ~y,
        yend = ~yend,
        type = "scatter",
        mode = "lines",
        text = ~hover_text,
        hoverinfo = "text",
        line = list(color = hit_col, width = 2, simplify = FALSE),
        showlegend = FALSE
      )
    }
  }

  p
}

aaSeqPanelPlotlySetFrameColors <- function(plot, frame_colors) {
  frame_colors <- frame_color_themes(frame_colors, FALSE)

  if (!is.null(plot$x$layout$shapes)) {
    for (i in seq_len(min(length(plot$x$layout$shapes), length(frame_colors)))) {
      plot$x$layout$shapes[[i]]$fillcolor <- frame_colors[i]
      plot$x$layout$shapes[[i]]$line$color <- frame_colors[i]
    }
  }

  if (length(plot$x$layoutAttrs) > 0 && !is.null(plot$x$layoutAttrs[[1]]$shapes)) {
    for (i in seq_len(min(length(plot$x$layoutAttrs[[1]]$shapes), length(frame_colors)))) {
      plot$x$layoutAttrs[[1]]$shapes[[i]]$fillcolor <- frame_colors[i]
      plot$x$layoutAttrs[[1]]$shapes[[i]]$line$color <- frame_colors[i]
    }
  }

  plot
}

aaSeqPanelPlotlySetRange <- function(plot, seq_length) {
  plot$x$data[[1]]$x <- c(1, seq_length)
  plot$x$data[[1]]$y <- c(-1, 2)

  if (length(plot$x$attrs) >= 1) {
    plot$x$attrs[[1]]$x <- c(1, seq_length)
    plot$x$attrs[[1]]$y <- c(-1, 2)
  }

  if (length(plot$x$layoutAttrs) == 0) {
    plot$x$layoutAttrs <- list(list())
  }
  if (is.null(plot$x$layoutAttrs[[1]]$xaxis)) plot$x$layoutAttrs[[1]]$xaxis <- list()
  if (is.null(plot$x$layoutAttrs[[1]]$yaxis)) plot$x$layoutAttrs[[1]]$yaxis <- list()
  plot$x$layoutAttrs[[1]]$xaxis$range <- c(1, seq_length)
  plot$x$layoutAttrs[[1]]$yaxis$range <- c(-1, 2)

  plot$x$layout$xaxis$range <- c(1, seq_length)
  plot$x$layout$yaxis$range <- c(-1, 2)

  if (!is.null(plot$x$layout$shapes)) {
    for (i in seq_along(plot$x$layout$shapes)) {
      plot$x$layout$shapes[[i]]$x1 <- seq_length
    }
  }

  plot
}

aaSeqPanelPlotlyTemplate <- function(frame_colors = "R") {
  frame_colors <- frame_color_themes(frame_colors, FALSE)
  shapes <- lapply(seq_len(3), function(i) {
    list(
      type = "rect",
      x0 = 1,
      x1 = 1,
      y0 = 2 - i,
      y1 = 3 - i,
      fillcolor = frame_colors[i],
      line = list(color = frame_colors[i], width = 0),
      layer = "below"
    )
  })

  p <- plotly::plotly_build(
    plotly::plot_ly(
      x = c(1, 1),
      y = c(-1, 2),
      type = "scatter",
      mode = "markers",
      hoverinfo = "none",
      showlegend = FALSE,
      marker = list(opacity = 0)
    ) %>%
      plotly::layout(
        margin = list(t = 0, r = 0, b = 40, l = 40, pad = 0),
        shapes = shapes,
        yaxis = list(
          autorange = FALSE,
          fixedrange = TRUE,
          range = c(-1, 2),
          zeroline = TRUE,
          showgrid = FALSE,
          showticklabels = FALSE,
          ticks = "",
          title = "frame"
        ),
        xaxis = list(
          autorange = FALSE,
          range = c(1, 1),
          showgrid = FALSE,
          title = list(text = "position [nt]")
        )
      )
  )

  aaSeqPanelPlotlySetRange(p, 1)
}

ntSeqPanelPlotlySetRange <- function(plot, seq_length) {
  plot$x$data[[1]]$x <- c(1, seq_length)
  plot$x$data[[1]]$y <- c(0, 1)

  if (length(plot$x$attrs) >= 1) {
    plot$x$attrs[[1]]$x <- c(1, seq_length)
    plot$x$attrs[[1]]$y <- c(0, 1)
  }

  if (length(plot$x$layoutAttrs) == 0) {
    plot$x$layoutAttrs <- list(list())
  }
  if (is.null(plot$x$layoutAttrs[[1]]$xaxis)) plot$x$layoutAttrs[[1]]$xaxis <- list()
  if (is.null(plot$x$layoutAttrs[[1]]$yaxis)) plot$x$layoutAttrs[[1]]$yaxis <- list()
  plot$x$layoutAttrs[[1]]$xaxis$range <- c(1, seq_length)
  plot$x$layoutAttrs[[1]]$yaxis$range <- c(0, 1)

  plot$x$layout$xaxis$range <- c(1, seq_length)
  plot$x$layout$yaxis$range <- c(0, 1)
  plot
}

ntSeqPanelPlotlyTemplate <- function() {
  p <- plotly::plotly_build(
    plotly::plot_ly(
      x = c(1, 1),
      y = c(0, 1),
      type = "scatter",
      mode = "markers",
      hoverinfo = "none",
      showlegend = FALSE,
      marker = list(opacity = 0)
    ) %>%
      plotly::layout(
        margin = list(t = 0, r = 0, b = 0, l = 0, pad = 0),
        yaxis = list(
          autorange = FALSE,
          fixedrange = TRUE,
          range = c(0, 1),
          showgrid = FALSE,
          showticklabels = FALSE,
          ticks = "",
          zeroline = FALSE,
          title = list(text = "")
        ),
        xaxis = list(
          autorange = FALSE,
          range = c(1, 1),
          showgrid = FALSE,
          showticklabels = FALSE,
          ticks = "",
          zeroline = FALSE,
          title = list(text = "")
        )
      )
  )
  ntSeqPanelPlotlySetRange(p, 1)
}

ntSeqPanelPlotly <- function(sequence, template = NULL) {
  seq_length <- nchar(sequence)

  if (inherits(template, "plotly")) {
    p <- template
    return(ntSeqPanelPlotlySetRange(p, seq_length))
  }

  p <- ntSeqPanelPlotlyTemplate()
  ntSeqPanelPlotlySetRange(p, seq_length)
}

seqPanelPlotTemplate <- function(theme_template) {
  suppressWarnings(ggplot(frame = "static")) +
    ylab("frame") +
    xlab("position [nt]") +
    theme_template +
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
  if (length(grl) <= 1L) {
    return(rep.int(1L, numExonsPerGroup(grl)))
  }
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

interval_overlaps_any <- function(query_start, query_end,
                                  subject_start, subject_end) {
  if (length(query_start) == 0L) return(logical())
  valid_query <- query_start <= query_end
  if (!any(valid_query) || length(subject_start) == 0L) {
    return(rep(FALSE, length(query_start)))
  }

  query_dt <- data.table(
    idx = seq_along(query_start),
    start = query_start,
    end = query_end
  )[valid_query]
  subject_dt <- unique(data.table(start = subject_start, end = subject_end))

  data.table::setkey(subject_dt, start, end)
  data.table::setkey(query_dt, start, end)

  hits <- data.table::foverlaps(query_dt, subject_dt, nomatch = 0L)
  out <- rep(FALSE, length(query_start))
  if (nrow(hits) > 0L) out[unique(hits$idx)] <- TRUE
  out
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
  overlap_type <- ifelse(viewMode == "tx" | collapse_intron_flank > 0, "within", "any")
  if (use_custom_region) {
    same_names <- names(custom_regions) %in% names(annotation)
    names(custom_regions)[same_names] <- paste(names(custom_regions)[same_names], "_1", sep="")
    overlaps_custom <- subsetByOverlaps(custom_regions, display_range,
                                        type = overlap_type)
  }

  display_names <- names(display_range)
  has_single_annotation <- length(annotation) == 1L &&
    length(display_range) == 1L &&
    identical(names(annotation), display_names)

  if (has_single_annotation) {
    annotation_shown <- annotation
  } else annotation_shown <- subsetByOverlaps(annotation, display_range,
                                      type = overlap_type)

  overlaps_tx <- NULL

  if (viewMode != "tx") {
    has_single_tx <- length(tx_annotation) == 1L &&
      length(display_range) == 1L &&
      identical(names(tx_annotation), display_names)

    if (has_single_tx) {
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

  return(gene_box_fix_overlaps(
    display_range,
    annotation_shown = annotation_shown,
    overlaps_tx = overlaps_tx,
    overlaps_custom = if (use_custom_region) overlaps_custom else NULL,
    custom_regions = custom_regions,
    collapse_intron_flank = collapse_intron_flank,
    frame_colors = frame_colors
  ))
}

flatten_gene_model_component <- function(grl, type) {
  if (is.null(grl) || length(grl) == 0) return(NULL)

  overlap_lengths <- width(grl@partitioning)
  grl_names <- grl@partitioning@NAMES
  if (is.null(grl_names)) grl_names <- as.character(seq_along(grl))
  names_grouping <- rep.int(grl_names, overlap_lengths)
  gr <- grl@unlistData
  gr@ranges@NAMES <- names_grouping
  gr$type <- rep.int(type, length(gr@ranges@start))
  gr
}

combine_gene_model_overlaps <- function(annotation_shown, overlaps_tx = NULL,
                                        overlaps_custom = NULL) {
  overlap_parts <- Filter(
    Negate(is.null),
    list(
      flatten_gene_model_component(annotation_shown, "cds"),
      flatten_gene_model_component(overlaps_tx, "utr"),
      flatten_gene_model_component(overlaps_custom, "cds")
    )
  )

  if (length(overlap_parts) == 0) {
    return(GenomicRanges::GRanges())
  }

  do.call(c, unname(overlap_parts))
}

gene_box_fix_overlaps <- function(display_range, annotation_shown, overlaps_tx,
                                  overlaps_custom, custom_regions,
                                  collapse_intron_flank, frame_colors) {
  overlaps <- combine_gene_model_overlaps(annotation_shown, overlaps_tx, overlaps_custom)
  if (length(overlaps) == 0) {
    result_dt <- data.table()
    lines_locations <- NULL
    return(list(result_dt, lines_locations))
  }
  overlaps$rel_frame_exon <- getRelativeFrames(overlaps)
  # overlaps <- subsetByOverlaps(overlaps, display_range)
  intersections <-  trimOverlaps(overlaps, display_range)
  intersections <- groupGRangesBy(intersections, paste0(names(intersections), "___", intersections$type))
  names(intersections) <- sub("___.*", "", names(intersections))

  locations_grl <- pmapToTranscriptF(intersections, display_range)
  intersection_lengths <- lengths(intersections)
  first_index <- cumsum(c(1L, head(intersection_lengths, -1L)))
  type_per_grl <- intersections@unlistData$type[first_index]
  locations_grl@unlistData$type <- rep.int(type_per_grl, times = lengths(locations_grl))
  locations <- unlistGrl(locations_grl)
  locations$rel_frame_orf <- getRelativeFrames(locations)

  cols <- colour_bars(locations, overlaps, display_range, frame_colors)
  locations_by_name <- groupGRangesBy(locations, locations@ranges@NAMES)
  group_layers_full <- geneTrackLayer(locations_by_name)
  group_lengths <- lengths(locations_by_name)
  group_first <- cumsum(c(1L, head(group_lengths, -1L)))
  layer_by_group <- group_layers_full[group_first]
  names(layer_by_group) <- locations_by_name@partitioning@NAMES
  layers <- unname(layer_by_group[names(locations)])
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

  location_ranges <- ranges(locations)
  rect_starts <- start(location_ranges)
  rect_ends <- end(location_ranges)
  gene_names <- names(location_ranges)

  blocks <- c(rect_starts, rect_ends)
  names(blocks) <- rep(gene_names, 2)
  blocks <- sort(blocks)
  lines_locations <- blocks[!(blocks %in% c(1, plot_width))]

  cap <- 0.00133*nchar(gene_names)
  labels_locations <- ((rect_starts + rect_ends) %/% 2L)
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
                          rect_starts = rect_starts,
                          rect_ends = rect_ends,
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
    seg_source <- unique(
      dt[no_ex > 1, .(gene_names, rect_starts, rect_ends, layers)]
    )
    data.table::setorder(seg_source, gene_names, rect_starts, rect_ends)

    did_collapse_introns <- FALSE
    if (collapse_intron_flank > 0 && nrow(seg_source) > 1L) {
      reduced_ranges <- GenomicRanges::reduce(
        GRanges(
          seg_source$gene_names,
          ranges = IRanges(seg_source$rect_starts, seg_source$rect_ends)
        )
      )
      did_collapse_introns <- length(reduced_ranges) != 1L
      if (did_collapse_introns) {
        seg_source <- data.table(
          gene_names = as.character(seqnames(reduced_ranges)),
          rect_starts = start(reduced_ranges),
          rect_ends = end(reduced_ranges)
        )
        unique_layers <- unique(dt[, .(gene_names, layers)])
        seg_source <- merge.data.table(
          seg_source,
          unique_layers,
          by = "gene_names",
          sort = FALSE
        )
        data.table::setorder(seg_source, gene_names, rect_starts, rect_ends)
      }
    }

    seg_dt <- seg_source[
      ,
      if (.N > 1L) {
        .(
          layers = layers[1L],
          rect_starts = head(rect_ends, -1L),
          rect_ends = tail(rect_starts, -1L),
          cols = "grey45",
          labels_locations = 0,
          hjusts = "center",
          type = "intron",
          no_ex = 1L
        )
      } else NULL,
      by = gene_names
    ]
    seg_dt <- seg_dt[rect_starts < rect_ends]
    if (did_collapse_introns && nrow(seg_dt) > 0L) {
      exon_intervals <- unique(dt[, .(rect_starts, rect_ends)])
      overlaps_other_exon <- interval_overlaps_any(
        seg_dt$rect_starts + 1L,
        seg_dt$rect_ends - 1L,
        exon_intervals$rect_starts,
        exon_intervals$rect_ends
      )
      seg_dt[
        (rect_ends - rect_starts > collapse_intron_flank * 2L) & !overlaps_other_exon,
        type := "intron_collapsed"
      ]
    }
    seg_dt <- seg_dt[, colnames(dt), with = FALSE]
    with_introns_dt <- rbindlist(list(dt, seg_dt))
    result_dt <- with_introns_dt
  }
  res_list <- list(result_dt, lines_locations)
  return(res_list)
}

move_utrs_touching_cds <- function(locations, type, plot_width) {
  utr_end_touch_cds <- type == "utr" & (end(locations)+1L) %in% start(locations)[type == "cds"]
  new_utr_ends <- end(locations)[utr_end_touch_cds] + 1L
  end(locations)[utr_end_touch_cds] <- pmin(plot_width, new_utr_ends)

  utr_start_touch_cds <- type == "utr" & (start(locations)-1L) %in% end(locations)[type == "cds"]
  new_utr_starts <- locations@ranges@start[utr_start_touch_cds] - 1L
  locations@ranges@start[utr_start_touch_cds] <- pmax(1L, new_utr_starts)
  return(locations)
}

alpha_normalize_colors <- function(colors, alpha = track_area_fill_alpha(), strength = 0.35) {
  color_matrix <- grDevices::col2rgb(colors)
  target <- round((color_matrix - (1 - alpha) * 255) / alpha)
  target[target < 0] <- 0
  target[target > 255] <- 255
  normalized <- round(color_matrix + strength * (target - color_matrix))
  apply(normalized, 2, function(v) {
    sprintf("#%02X%02X%02X", v[[1]], v[[2]], v[[3]])
  })
}

frame_color_themes <- function(theme, with_alpha = FALSE) {
  if (theme == "R") {
    colors <- rc_rgb(FALSE)
  } else if (theme == "Color_blind") {
    colors <- c("#5056A0", "#F9AA2A", "#AED5EB")
  } else stop("theme must be either R or Color blind")

  if (with_alpha) colors <- alpha_normalize_colors(colors)
  colors
}

rc_rgb <- function(with_alpha = TRUE) {
  rgb_colors <- c("#F8766D", "#00BA38", "#619CFF")
  if (with_alpha) rgb_colors <- alpha_normalize_colors(rgb_colors)
  rgb_colors
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
geneModelPanelPlotlyTemplate <- function() {
  plotly::plot_ly(
    x = numeric(0),
    y = numeric(0),
    type = "scatter",
    mode = "markers",
    hoverinfo = "skip",
    showlegend = FALSE
  ) %>%
    plotly::layout(
      showlegend = FALSE,
      margin = list(l = 0, r = 0, b = 0, t = 0),
      paper_bgcolor = "white",
      plot_bgcolor  = "white",
      hovermode = "closest"
    )
}

geneModelPanelPlotly <- function(dt, template = geneModelPanelPlotlyTemplate()) {
  # If no annotation given, return blank plotly
  if (is.null(dt) || nrow(dt) == 0) {
    p <- template
    p$x$data <- list(list(
      x = 0,
      y = 0,
      type = "scatter",
      mode = "markers",
      marker = list(opacity = 0),
      hoverinfo = "skip",
      showlegend = FALSE
    ))
    return(p %>%
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

  all_shapes <- list()
  p <- template %>%
    plotly::layout(
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
      )
    )

  # --- Introns as line shapes rendered behind boxes ---
  if (nrow(seg_dt) > 0) {
    seg_shapes <- lapply(seq_len(nrow(seg_dt)), function(i) {
      list(
        type = "line",
        layer = "below",
        xref = "x", yref = "y",
        x0 = seg_dt$rect_starts[i],
        x1 = seg_dt$rect_ends[i],
        y0 = 0.5 - seg_dt$layers[i],
        y1 = 0.5 - seg_dt$layers[i],
        line = list(color = "grey45", width = 2),
        opacity = 0.6
      )
    })
    all_shapes <- c(all_shapes, seg_shapes)

    p <- p %>%
      plotly::add_markers(
        data = seg_dt,
        x = ~(rect_starts + rect_ends) / 2,
        y = ~(0.5 - layers),
        marker = list(opacity = 0, size = 10),
        text = ~gene_names,
        hovertext = ~gene_names,
        hoverinfo = "text",
        showlegend = FALSE
      )

    # Collapsed intron label "..."
    intron_flank_coords <- seg_dt[type %in% c("intron_collapsed")]
    if (nrow(intron_flank_coords) > 0) {
      intron_flank_coords[, `:=`(
        xmid = (rect_starts + rect_ends) / 2,
        ymid = 0.5 - layers,
        tooltip = "Collapsed Intron"
      )]
      p <- p %>%
        plotly::add_text(
          data = intron_flank_coords,
          x = ~xmid, y = ~ymid,
          text = "...",
          textposition = "middle center",
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
    all_shapes <- c(all_shapes, rect_shapes)

    # Add an invisible scatter trace to provide per-rectangle hover tooltips
    box_dt[, `:=`(
      xmid = (rect_starts + rect_ends) / 2,
      ymid = (ymin + ymax) / 2
    )]

    p <- p %>%
      plotly::add_markers(
        data = box_dt,
        x = ~xmid, y = ~ymid,
        marker = list(opacity = 0, size = 12),
        text = ~gene_names,
        hovertext = ~gene_names,
        hoverinfo = "text"
      )

    # Gene name labels per (layer, gene_names) at mean(labels_locations)
    lab_dt <- box_dt[, .(
      layers = layers[1],
      labels_locations = mean(labels_locations)
    ), by = gene_names]

    p <- p %>%
      plotly::add_text(
        data = lab_dt,
        x = ~labels_locations,
        y = ~(0.50 - layers),
        text = ~gene_names,
        textposition = "middle center",
        textfont = list(color = "black", size = 16),
        hoverinfo = "skip",
        name = "id",
        legendgroup = "id",
        showlegend = TRUE
      )
  }
  if (length(all_shapes) > 0) {
    p <- p %>% plotly::layout(shapes = all_shapes)
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
