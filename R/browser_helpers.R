browser_tx_seqs_getSeq <- function(reference_sequence, display_range,
                                   keep.names = FALSE) {
  if (is(display_range, "GRanges")) {
    display_range <- GenomicRanges::GRangesList(tx = display_range)
  }
  stopifnot(is(display_range, "GRangesList"))

  tx_count <- length(display_range)
  if (tx_count == 0L) return(Biostrings::DNAStringSet())

  tx_names <- names(display_range)
  if (is.null(tx_names)) tx_names <- as.character(seq_len(tx_count))

  if (tx_count == 1L) {
    gr <- display_range[[1]]
    tx_strand <- as.character(GenomicRanges::strand(gr)[1])

    if (length(gr) > 1L && identical(tx_strand, "-")) {
      starts <- BiocGenerics::start(gr)
      if (!all(starts[-1L] >= starts[-length(starts)])) {
        gr <- gr[order(starts, BiocGenerics::end(gr))]
      }
    }

    GenomicRanges::strand(gr) <- "+"
    exon_seqs <- if (is(reference_sequence, "FaFile")) {
      Rsamtools::scanFa(reference_sequence, param = gr)
    } else {
      Biostrings::getSeq(reference_sequence, gr)
    }
    tx_string <- if (length(exon_seqs) == 1L) {
      exon_seqs[[1]]
    } else {
      Biostrings::DNAString(paste0(as.character(exon_seqs), collapse = ""))
    }

    if (identical(tx_strand, "-")) {
      tx_string <- Biostrings::reverseComplement(tx_string)
    }
    return(tx_string)
  }

  exon_counts <- lengths(display_range)

  gr <- unlist(display_range, use.names = FALSE)
  group_starts <- cumsum(c(1L, head(exon_counts, -1L)))
  tx_strands <- as.character(GenomicRanges::strand(gr)[group_starts])
  tx_ends <- cumsum(exon_counts)

  if (tx_count == 1L && exon_counts[[1]] > 1L && identical(tx_strands[[1]], "-")) {
    ord <- order(BiocGenerics::start(gr), BiocGenerics::end(gr))
    gr <- gr[ord]
  } else if (any(exon_counts > 1L & tx_strands == "-")) {
    exon_idx <- seq_along(gr)
    tx_starts <- c(1L, head(tx_ends, -1L) + 1L)
    for (i in which(exon_counts > 1L & tx_strands == "-")) {
      grp_idx <- tx_starts[[i]]:tx_ends[[i]]
      exon_idx[grp_idx] <- grp_idx[order(
        BiocGenerics::start(gr)[grp_idx],
        BiocGenerics::end(gr)[grp_idx]
      )]
    }
    gr <- gr[exon_idx]
  }

  # Fetch genomic sequence directly and apply transcript-level reverse
  # complement after exon concatenation for minus-strand transcripts.
  GenomicRanges::strand(gr) <- "+"
  exon_seqs <- if (is(reference_sequence, "FaFile")) {
    Rsamtools::scanFa(reference_sequence, param = gr)
  } else {
    Biostrings::getSeq(reference_sequence, gr)
  }

  if (all(exon_counts == 1L)) {
    tx_strings <- as.character(exon_seqs)
  } else {
    tx_starts <- c(1L, head(tx_ends, -1L) + 1L)
    exon_strings <- as.character(exon_seqs)
    tx_strings <- character(tx_count)
    for (i in seq_len(tx_count)) {
      tx_strings[[i]] <- paste0(exon_strings[tx_starts[[i]]:tx_ends[[i]]], collapse = "")
    }
  }

  minus_tx <- tx_strands == "-"
  if (any(minus_tx)) {
    tx_strings[minus_tx] <- vapply(
      tx_strings[minus_tx],
      function(x) as.character(Biostrings::reverseComplement(Biostrings::DNAString(x))),
      character(1)
    )
  }

  tx_seqs <- Biostrings::DNAStringSet(tx_strings)

  if (isTRUE(keep.names)) {
    tx_seqs@ranges@NAMES <- tx_names
  } else {
    tx_seqs@ranges@NAMES <- NULL
  }
  tx_seqs
}

multiOmicsPlot_bottom_panels <- function(reference_sequence, display_range, annotation,
                                         start_codons, stop_codons, custom_motif,
                                         custom_regions, viewMode,
                                         tx_annotation = NULL, collapse_intron_flank = 100,
                                         frame_colors = "R",
                                         templates = NULL) {
  force(display_range)
  # Get sequence and create basic seq panel
  target_seq <- browser_tx_seqs_getSeq(reference_sequence, display_range)
  seq_panel_hits <- createSeqPanelPattern(target_seq, start_codons = start_codons,
                                          stop_codons = stop_codons, custom_motif = custom_motif)
  seq_aa_panel <- plotAASeqPanelPlotly(
    seq_panel_hits,
    target_seq,
    frame_colors = frame_colors,
    template = templates$aa_seq_panel_plotly
  )
  # Get the panel for the annotation track
  gene_model_panel_dt <- createGeneModelPanel(display_range, annotation,
                                           tx_annotation = tx_annotation,
                                           custom_regions = custom_regions,
                                           viewMode = viewMode, collapse_intron_flank,
                                           frame_colors = frame_colors)
  lines <- gene_model_panel_dt[[2]]
  layers <- if (nrow(gene_model_panel_dt[[1]]) == 0) 1L else max(gene_model_panel_dt[[1]]$layers)

  gene_model_panel <- geneModelPanelPlotly(gene_model_panel_dt[[1]])
  seq_nt_panel <- ntSeqPanelPlotly(
    target_seq,
    template = templates$nt_seq_panel_plotly
  )
  return(list(seq_panel = seq_aa_panel, seq_nt_panel = seq_nt_panel,
              gene_model_panel_dt = gene_model_panel_dt[[1]],
              gene_model_panel = gene_model_panel, frame_colors = frame_colors,
              lines = lines, target_seq = target_seq, annotation_layers = layers))
}

multiOmicsPlot_all_track_plots <- function(profiles, withFrames, frame_colors, colors,
                                           ylabels, ylabels_full_name, lines,
                                           frames_type, summary_track, summary_track_type,
                                           BPPARAM, templates = NULL) {
  force(colors)
  force(lines)
  force(ylabels)

  if (frames_type == "animate") {
    plots <- list(getPlotAnimate(rbindlist(profiles, idcol = "file"), withFrames = withFrames[1],
                                colors = colors[1], frame_colors = frame_colors,
                                ylabels = ylabels[1], lines = lines))
  } else {
    plots <- build_browser_track_plots(
      profiles,
      withFrames,
      frame_colors,
      colors,
      ylabels,
      ylabels_full_name,
      lines,
      frames_type,
      BPPARAM,
      templates = templates
    )
  }

  nplots <- ifelse(frames_type == "animate", 1, length(plots))
  if (summary_track) {
    nplots <- nplots + 1
    plots <- make_summary_track(profiles, plots, withFrames, frame_colors, colors,
                                lines, summary_track_type, nplots, templates = templates)
  }
  return(list(plots = plots, nplots = nplots, track_type = frames_type))
}

multiOmicsPlot_all_profiles <- function(display_range, reads, kmers,
                                        kmers_type, frames_type, frames_subset,
                                        withFrames, log_scale, BPPARAM) {
  force(display_range)
  force(reads)
  force(kmers)
  force(kmers_type)
  force(frames_type)
  force(withFrames)
  force(log_scale)
  force(frames_subset)
  profiles <- build_browser_profiles(
    display_range,
    reads,
    kmers,
    kmers_type,
    frames_type,
    frames_subset,
    withFrames,
    log_scale,
    BPPARAM
  )
  return(profiles)
}

multiOmicsPlot_complete_plot <- function(track_panel, bottom_panel, display_range,
                                         proportions, seq_render_dist,
                                         display_sequence,
                                         aa_letter_code, input_id, plot_name,
                                         plot_title,  width, height, export.format,
                                         zoom_range = NULL, frame_colors = "R") {
  print("Merging coverage tracks and bottom")
  track_plots <- browser_plots_highlighted(track_panel$plots, zoom_range)
  bottom_plots <- bottom_panel$bottom_plots
  plots <- c(track_plots, bottom_plots)

  old <- TRUE
  if (old) {
    multiomics_plot <- suppressWarnings(subplot(plots,
                                                margin = 0,
                                                nrows = length(plots),
                                                heights = proportions,
                                                shareX = TRUE,
                                                titleY = TRUE, titleX = TRUE))
  } else {
    multiomics_plot <- suppressWarnings(fast_subplot_shared_x(plots,
                                                              margin = 0,
                                                              nrows = length(plots),
                                                              heights = proportions,
                                                              shareX = TRUE,
                                                              titleY = TRUE,
                                                              titleX = TRUE))
  }



  if (isTruthy(display_sequence)) {
    nt_seq_y_index <- length(plots) - bottom_panel$ncustom - 3
    multiomics_plot <- addJSrender(multiomics_plot, bottom_panel$target_seq,
                                   nt_seq_y_index, seq_render_dist,
                                   aa_letter_code, input_id, bottom_panel$frame_colors)
  }

  return(browser_plot_final_layout_polish(multiomics_plot, plot_name, display_range,
                                          width, height, export.format, plot_title,
                                          zoom_range, proportions,
                                          apply_line_desimplify = identical(track_panel$track_type, "animate")))
}

#' @importFrom Biostrings nchar translate
#'
multiOmicsPlot_internal <- function(display_range, df, annotation = "cds", reference_sequence = findFa(df),
                                    reads = outputLibs(df, type = "pshifted", output.mode = "envirlist",
                                                       naming = "full", BPPARAM = BiocParallel::SerialParam()),
                                    viewMode = c("tx", "genomic")[1],
                                    custom_regions = NULL,
                                    leader_extension = 0, trailer_extension = 0,
                                    withFrames = libraryTypes(df, uniqueTypes = FALSE) %in% c("RFP", "RPF", "LSU", "TI"),
                                    frames_type = "lines", colors = NULL, kmers = NULL, kmers_type = c("mean", "sum")[1],
                                    ylabels = bamVarName(df), lib_to_annotation_proportions = c(0.8,0.2),lib_proportions = NULL,
                                    annotation_proportions = NULL,width = NULL, height = NULL,
                                    plot_name = "default", plot_title = NULL,
                                    display_sequence = c("both","nt", "aa", "none")[1], seq_render_dist = 100,
                                    aa_letter_code = c("one_letter", "three_letters")[1],
                                    annotation_names = NULL, start_codons = "ATG", stop_codons = c("TAA", "TAG", "TGA"),
                                    custom_motif = NULL, log_scale = FALSE, BPPARAM = BiocParallel::SerialParam(),
                                    input_id = "", summary_track = FALSE,
                                    summary_track_type = frames_type, export.format = "svg", frames_subset = "all",
                                    zoom_range = NULL, tx_annotation = NULL, collapse_intron_flank = 100,
                                    frame_colors = "R", phyloP = FALSE, mapability = FALSE) {

  multiOmicsController()
  # Get Bottom annotation and sequence panels
  bottom_panel <- multiOmicsPlot_bottom_panels(reference_sequence, display_range, annotation,
                                               start_codons, stop_codons, custom_motif,
                                               custom_regions, viewMode,
                                               tx_annotation,
                                               collapse_intron_flank,
                                               frame_colors)
  custom_bigwig_panels <- custom_seq_track_panels(df, display_range, phyloP, mapability)
  bottom_panel <- c(bottom_panel, annotation_list, custom_bigwig_panels,
                    ncustom = length(custom_bigwig_panels[[1]]))

  bottom_panel$bottom_plots <- bottom_plots_to_plotly(bottom_panel)


  # Get NGS data track panels
  multiOmicsControllerView()
  profiles <- multiOmicsPlot_all_profiles(display_range, reads, kmers,
                                          kmers_type, frames_type, frames_subset,
                                          withFrames, log_scale, BPPARAM)

  track_panel <- multiOmicsPlot_all_track_plots(profiles, withFrames,
                                          frame_colors, colors, ylabels,
                                          ylabels_full_name, bottom_panel$lines,
                                          frames_type,
                                          summary_track, summary_track_type,
                                          BPPARAM)
  return(multiOmicsPlot_complete_plot(track_panel, bottom_panel, display_range,
                                      proportions, seq_render_dist,
                                      display_sequence,
                                      aa_letter_code, input_id, plot_name,
                                      plot_title,  width, height, export.format,
                                      zoom_range, frame_colors))
}

genomic_string_to_grl <- function(genomic_string, display_region, max_size = 1e6,
                                  viewMode, extendLeaders, extendTrailers,
                                  collapsed_introns_width = 0, type = "region") {
  input_given <- !is.null(genomic_string) && genomic_string != ""
  if (input_given) {
    gr <- try(as(unlist(strsplit(genomic_string, ";")), "GRanges"))
    if (is(gr, "GRanges")) {
      if (any(start(gr) < 1)) stop("Position 1 is minimum position to show on a chromosome! (input ",
                                   paste(start(gr), collapse = ", "), ")")

      suppressWarnings(try(seqlevelsStyle(gr) <- seqlevelsStyle(display_region)[1], silent = TRUE))

      if (any(!(as.character(seqnames(gr)) %in% seqnames(seqinfo(display_region)))))
        stop("Invalid chromosome selected!")
      display_region <- GRangesList(Region = gr)
    } else stop("Malformed genomic ", type, ": format: chr1:39517672-39523668:+;chr1:39527673-39520669:+")
  }
  extension_size <- extendLeaders + extendTrailers
  true_sized_grl <- if (!viewMode) {
    display_region
  }  else {
    if (collapsed_introns_width > 0) {
      print("Collapsing introns")
      display_region <- exonsWithPseudoIntronsPerGroup(display_region, collapsed_introns_width)
    } else display_region <- flankPerGroup(display_region)
  }

  size <- widthPerGroup(true_sized_grl, FALSE)
  if (size > max_size) stop("Only up to ", round(max_size/1e6, 3) ," million bases can be shown, input: ",
                            round(size/1e6, 3), " million bases")
  return(display_region)
}

get_zoom_range <- function(zoom_range, display_region, max_size,
                           viewMode, leader_extension, trailer_extension,
                           zoom_range_flank = 10) {
  if (!is.null(zoom_range)) {
    valid_zoom <- is.character(zoom_range) && nchar(zoom_range) > 0
    if (valid_zoom) {
      count_colon <- stringr::str_count(zoom_range, ":")
      tx_coord_interval <- count_colon == 1
      geomic_coord_interval <- count_colon > 1
      if (tx_coord_interval) {
        zoom_interval <- as.numeric(unlist(strsplit(zoom_range, ":")))
        stopifnot(length(zoom_interval) == 2 && zoom_interval[1] <= zoom_interval[2])
        zoom_range <- zoom_interval
      } else if (geomic_coord_interval) {
        gr <- genomic_string_to_grl(zoom_range, display_region, max_size,
                                    viewMode, leader_extension, trailer_extension,
                                    "zoom region")
        display_range_zoom <- display_region
        if (viewMode) {
          display_range_zoom <- flankPerGroup(display_range_zoom)
          gr <- flankPerGroup(gr)
        }
        if (!is.null(leader_extension) && is.numeric(leader_extension) && leader_extension != 0)
          display_range_zoom <- extendLeaders(display_range_zoom, leader_extension)
        if (!is.null(trailer_extension) && is.numeric(trailer_extension) &&  trailer_extension != 0)
          display_range_zoom <- extendTrailers(display_range_zoom, trailer_extension)
        ir <- suppressWarnings(pmapToTranscriptF(gr, display_range_zoom))
        if (as.numeric(width(ir)) > 0) {
          zoom_range <- c(max(as.numeric(start(ir))[1] - zoom_range_flank, 1),
                          min(as.numeric(end(ir))[1] + zoom_range_flank,
                              widthPerGroup(display_range_zoom, FALSE)))
        } else {
          zoom_range <- numeric(0)
          if (!viewMode) {
            attr(zoom_range, "message") <-
            "Zoom range not overlapping displayed region, did you mean to use genomic coordinates?"
          } else {
            attr(zoom_range, "message") <-
              "Zoom range not overlapping displayed region."
          }
          warning(attr(zoom_range, "message"))
        }
      }
    }
  }
  if (!is.numeric(zoom_range)) zoom_range <- numeric(0)
  return(zoom_range)
}

browser_plots_highlighted <- function(plots, zoom_range, color = "rgba(255, 255, 102, 0.18)") {
  if (!is.null(zoom_range) && length(zoom_range) == 2) {
    plots_highlighted <- lapply(plots, function (p) {
      p$x$layout$shapes <- c(
        p$x$layout$shapes,
        list(list(
          type = "rect",
          yref = "paper",  # "paper" ensures full Y-axis coverage
          x0 = zoom_range[1], x1 = zoom_range[2],
          y0 = 0, y1 = 1,  # Full height
          fillcolor = color,  # Light yellow
          line = list(width = 0) # Remove border
        ))
      )
      return(p)
    })
    plots <- plots_highlighted
  }
  return(plots)
}

browser_input_or_default <- function(input, name, default = NULL) {
  value <- tryCatch(input[[name]], error = function(e) NULL)
  if (is.null(value) || length(value) == 0) default else value
}

browser_legend_cleanup <- function(plot) {
  `%||%` <- function(x, y) if (is.null(x)) y else x
  seen_legends <- character()

  for (i in seq_along(plot$x$data)) {
    trace <- plot$x$data[[i]]
    legend_key <- trace$legendgroup %||% trace$name

    if (isTRUE(trace$showlegend) && !is.null(legend_key) && nzchar(legend_key)) {
      if (legend_key %in% seen_legends) {
        plot$x$data[[i]]$showlegend <- FALSE
      } else {
        seen_legends <- c(seen_legends, legend_key)
      }
    }
  }

  plot
}

observatory_selection_cache_key <- function(library_selections,
                                            library_selection_labels = NULL) {
  if (is.null(library_selections) || length(library_selections) == 0) return("")

  selection_ids <- names(library_selections)
  if (is.null(selection_ids)) selection_ids <- as.character(seq_along(library_selections))

  parts <- vapply(seq_along(library_selections), function(i) {
    selection_id <- selection_ids[[i]]
    label <- ""
    if (!is.null(library_selection_labels) && !is.null(library_selection_labels[[selection_id]])) {
      label <- as.character(library_selection_labels[[selection_id]])
    }
    runs <- as.character(library_selections[[i]])
    runs <- sort(unique(runs[!is.na(runs) & nzchar(runs)]))
    paste(selection_id, label, paste(runs, collapse = ","), sep = ":")
  }, character(1))

  paste(parts, collapse = "|group|")
}

hash_strings_browser <- function(input, dff, ciw = input$collapsed_introns_width) {
  full_names <- runIDs(dff)
  if (all(full_names == "")) full_names <- orfik_name_decider(dff, naming = "full")

  hash_bottom <- paste(browser_input_or_default(input, "tx", ""),
                       browser_input_or_default(input, "other_tx", FALSE),
                       browser_input_or_default(input, "add_uorfs", FALSE),
                       browser_input_or_default(input, "add_translon", FALSE),
                       browser_input_or_default(input, "add_translons_transcode", FALSE),
                       browser_input_or_default(input, "extendTrailers", 0),
                       browser_input_or_default(input, "extendLeaders", 0),
                       browser_input_or_default(input, "genomic_region", ""),
                       browser_input_or_default(input, "viewMode", FALSE),
                       ciw,
                       browser_input_or_default(input, "colors", ""),
                       browser_input_or_default(input, "customSequence", ""),
                       browser_input_or_default(input, "phyloP", FALSE),
                       browser_input_or_default(input, "mapability", FALSE),
                       collapse = "|_|")
  # Until plot and coverage is split (bottom must be part of browser hash)
  hash_browser <- paste(hash_bottom,
                        full_names,
                        browser_input_or_default(input, "plot_export_format", "svg"),
                        browser_input_or_default(input, "summary_track", FALSE),
                        browser_input_or_default(input, "summary_track_type", "lines"),
                        browser_input_or_default(input, "kmer", 1),
                        browser_input_or_default(input, "frames_type", "lines"),
                        browser_input_or_default(input, "withFrames", TRUE),
                        browser_input_or_default(input, "log_scale", FALSE),
                        browser_input_or_default(input, "zoom_range", ""),
                        browser_input_or_default(input, "frames_subset", "all"),
                        browser_input_or_default(input, "unique_align", FALSE),
                        collapse = "|_|")
  hash_expression <- paste(full_names,
                           browser_input_or_default(input, "tx", ""),
                           browser_input_or_default(input, "expression_plot", FALSE),
                           browser_input_or_default(input, "extendTrailers", 0),
                           browser_input_or_default(input, "extendLeaders", 0),
                           collapse = "|_|")
  hash_strings <- list(hash_bottom = hash_bottom, hash_browser = hash_browser,
                       hash_expression = hash_expression)
  stopifnot(all(lengths(hash_strings) == 1))
  return(hash_strings)
}

# Lock yaxis domains from a proportions vector (top -> bottom)
lock_yaxis_domains_by_proportions <- function(p, proportions, gap = 0, fixed = TRUE, layer_below = TRUE,
                                              uirevision = TRUE, margins = list(t = 30, b = 60, l = 60, r = 10)) {
  stopifnot(length(proportions) >= 1)
  n <- length(proportions)

  # normalize proportions and compute absolute heights after gaps
  prop <- proportions / sum(proportions)
  total_gap <- (n - 1) * gap
  avail <- 1 - total_gap
  heights <- prop * avail

  # build yaxisN domains (Plotly expects y=0 bottom, y=1 top; yaxis is the TOP row)
  layout_updates <- list()
  current_top <- 1.0
  for (i in seq_len(n)) {
    h  <- heights[i]
    y1 <- current_top
    y0 <- y1 - h
    ax <- if (i == 1) "yaxis" else paste0("yaxis", i)

    ya <- list(domain = c(y0, y1))
    if (fixed)        ya$fixedrange <- FALSE
    if (layer_below)  ya$layer <- "below traces"

    layout_updates[[ax]] <- ya

    # move down for next row (add gap below this row)
    current_top <- y0 - gap
  }

  # apply with do.call (no tidy-eval)
  p <- do.call(
    plotly::layout,
    c(
      list(p, uirevision = uirevision, margin = margins),
      layout_updates
    )
  )
#   onRender(p, "
# function(el, x){
#   const gd = document.getElementById(el.id);
#   let guarding = false;
#   gd.on('plotly_relayout', ev => {
#     if (guarding) return;
#     // if the y range changed, clamp the lower bound to 0
#     const y0 = ev['yaxis.range[0]'];
#     const y1 = ev['yaxis.range[1]'];
#     const y = ev['yaxis.range'];
#     if (y0 !== undefined || y1 !== undefined || Array.isArray(y)) {
#       const newMax = Array.isArray(y) ? y[1] : (y1 !== undefined ? y1 : gd.layout.yaxis.range[1]);
#       guarding = true;
#       Plotly.relayout(gd, {'yaxis.range': [0, newMax]}).then(() => guarding = false);
#     }
#   });
# }
# ")
  p
}

browser_plot_final_layout_polish <- function(multiomics_plot,
                                             plot_name,
                                             display_range,
                                             width,
                                             height,
                                             export.format,
                                             plot_title,
                                             zoom_range,
                                             proportions,
                                             apply_line_desimplify = FALSE) {
  multiomics_plot <- remove_y_axis_zero_tick_js(multiomics_plot)

  xaxis_names <- names(multiomics_plot$x$layout)[grepl("^xaxis[0-9]*$", names(multiomics_plot$x$layout))]
  bottom_xaxis <- "xaxis"
  if (length(xaxis_names) > 0) {
    `%||%` <- function(x, y) if (is.null(x)) y else x

    axis_bottom <- function(axis_name) {
      axis <- multiomics_plot$x$layout[[axis_name]]
      anchor <- axis$anchor %||% "y"
      yaxis_name <- if (identical(anchor, "y")) "yaxis" else paste0("yaxis", sub("^y", "", anchor))
      yaxis <- multiomics_plot$x$layout[[yaxis_name]]
      domain <- yaxis$domain %||% c(0, 1)
      domain[[1]]
    }

    bottom_xaxis <- xaxis_names[[which.min(vapply(xaxis_names, axis_bottom, numeric(1)))]]

    for (axis_name in xaxis_names) {
      axis <- multiomics_plot$x$layout[[axis_name]]
      if (is.null(axis)) axis <- list()
      axis$visible <- identical(axis_name, bottom_xaxis)
      axis$showticklabels <- identical(axis_name, bottom_xaxis)
      axis$tickfont <- list(size = 16)
      axis$title <- if (identical(axis_name, bottom_xaxis)) {
        list(text = "position [nt]", font = list(size = 22))
      } else {
        list(text = "")
      }
      multiomics_plot$x$layout[[axis_name]] <- axis
    }
  }

  multiomics_plot <- multiomics_plot %>%
    plotly::layout(
      legend = list(
        x = 1.02,
        xanchor = "left",
        y = 0.92,
        yanchor = "top",
        orientation = "v"
      )
    )
  multiomics_plot$x$layout$legend <- list(
    x = 1.02,
    xanchor = "left",
    y = 0.92,
    yanchor = "top",
    orientation = "v"
  )

  filename <- ifelse(plot_name == "default", names(display_range), plot_name)
  multiomics_plot <- addToImageButtonOptions(multiomics_plot, filename,
                                             width, height, format = export.format)
  if (!is.null(plot_title)) multiomics_plot <- multiomics_plot %>%
    plotly::layout(title = plot_title)
  # Lock proportions on zoom out
  # multiomics_plot <- lock_yaxis_domains_by_proportions(multiomics_plot, proportions, gap = 0)
  if (!is.null(zoom_range) && length(zoom_range) == 2) {
    # Apply the initial zoom to every shared x-axis so added tracks do not
    # reset the visible axis back to the full range.
    for (axis_name in unique(c(bottom_xaxis, xaxis_names))) {
      axis <- multiomics_plot$x$layout[[axis_name]]
      if (is.null(axis)) axis <- list()
      axis$range <- zoom_range
      multiomics_plot$x$layout[[axis_name]] <- axis
    }
  }
  max_x <- as.numeric(widthPerGroup(display_range, FALSE))
  if (length(max_x) > 1) max_x <- max_x[[1]]
  multiomics_plot <- addBrowserXRangeClamp(multiomics_plot, min_x = 1, max_x = max_x)
  multiomics_plot <- browser_legend_cleanup(multiomics_plot)
  multiomics_plot <- addColumnsZoomSwitch(multiomics_plot)
  if (isTRUE(apply_line_desimplify)) {
    return(lineDeSimplify(multiomics_plot))
  }
  multiomics_plot
}

fast_subplot_shared_x <- function(...,
                                  nrows = NULL,
                                  heights = NULL,
                                  shareX = FALSE,
                                  titleY = TRUE,
                                  titleX = TRUE,
                                  margin = 0) {
  `%||%` <- function(x, y) if (is.null(x)) y else x

  args <- list(...)
  plots <- if (length(args) == 1L && is.list(args[[1]]) && !
               inherits(args[[1]], "plotly")) {
    args[[1]]
  } else {
    args
  }

  if (!length(plots)) stop("No plots supplied.")
  if (!all(vapply(plots, inherits, logical(1), what = "plotly"))) {
    stop("All inputs must be plotly objects.")
  }
  if (!isTRUE(shareX)) {
    stop("fast_subplot_shared_x currently only supports shareX = TRUE.")
  }

  n <- length(plots)
  if (is.null(nrows)) nrows <- n
  if (length(nrows) != 1L || is.na(nrows) || nrows != n) {
    stop("fast_subplot_shared_x currently only supports one plot per row.")
  }

  if (is.null(heights)) heights <- rep(1, n)
  if (length(heights) != n) stop("`heights` must have one value per plot.")
  if (n > 1L && (length(margin) != 1L || !is.numeric(margin) || is.na(margin)
                 || margin < 0)) {
    stop("`margin` must be a single non-negative number.")
  }

  built <- lapply(plots, function(p) {
    is_materialized <- !is.null(p$x$data) &&
      !is.null(p$x$layout) &&
      length(p$x$layoutAttrs %||% list()) == 0L
    if (is_materialized) {
      p
    } else {
      suppressWarnings(plotly::plotly_build(p))
    }
  })
  layouts <- lapply(built, function(p) p$x$layout %||% list())

  total_gap <- margin * max(0, n - 1L)
  if (total_gap >= 1) stop("`margin` is too large for the number of rows.")
  heights <- heights / sum(heights) * (1 - total_gap)

  domains <- vector("list", n)
  top <- 1
  for (i in seq_len(n)) {
    bottom <- top - heights[[i]]
    domains[[i]] <- c(bottom, top)
    top <- bottom - margin
  }

  axis_ref <- function(prefix, i) {
    if (i == 1L) prefix else paste0(prefix, i)
  }

  remap_axis_ref <- function(ref, target, axis_prefix) {
    if (is.null(ref)) return(target)

    ref_chr <- as.character(ref)
    if (identical(ref_chr, "paper")) return("paper")
    if (grepl(paste0("^", axis_prefix, "[0-9]* domain$"), ref_chr)) {
      return(paste0(target, " domain"))
    }
    if (grepl(paste0("^", axis_prefix, "[0-9]*$"), ref_chr)) {
      return(target)
    }
    ref_chr
  }

  out <- built[[1]]
  trace_counts <- vapply(built, function(p) length(p$x$data %||% list()), integer(1))
  attrs_counts <- vapply(built, function(p) length(p$x$attrs %||% list()), integer(1))
  shape_counts <- vapply(layouts, function(lay) length(lay$shapes %||% list()), integer(1))
  annotation_counts <- vapply(layouts, function(lay) length(lay$annotations %||% list()), integer(1))
  image_counts <- vapply(layouts, function(lay) length(lay$images %||% list()), integer(1))

  out$x$data <- vector("list", sum(trace_counts))
  out$x$attrs <- vector("list", sum(attrs_counts))
  out$x$layout <- out$x$layout %||% list()
  out$x$layoutAttrs <- list()
  out$x$layout$shapes <- vector("list", sum(shape_counts))
  out$x$layout$annotations <- vector("list", sum(annotation_counts))
  out$x$layout$images <- vector("list", sum(image_counts))
  out$x$visdat <- NULL
  out$x$cur_data <- NULL

  shared_xaxis <- layouts[[n]]$xaxis %||% layouts[[1]]$xaxis %||% list()
  bottom_yaxis_name <- paste0("yaxis", if (n == 1L) "" else n)
  shared_xaxis$domain <- c(0, 1)
  shared_xaxis$anchor <- sub("^yaxis", "y", bottom_yaxis_name)
  if (!isTRUE(titleX)) shared_xaxis$title <- list(text = "")
  out$x$layout$xaxis <- shared_xaxis

  out$x$layout$margin <- layouts[[1]]$margin %||% list(l = 0, r = 0, t = 0, b = 0, pad = 0)

  trace_index <- 1L
  attrs_index <- 1L
  shape_index <- 1L
  annotation_index <- 1L
  image_index <- 1L

  for (i in seq_len(n)) {
    p <- built[[i]]
    layout_i <- layouts[[i]]

    yaxis_name <- axis_ref("yaxis", i)
    ytrace_name <- axis_ref("y", i)

    yaxis <- layout_i$yaxis %||% list()
    yaxis$domain <- domains[[i]]
    yaxis$anchor <- "x"
    if (!isTRUE(titleY)) yaxis$title <- list(text = "")
    out$x$layout[[yaxis_name]] <- yaxis

    traces_i <- p$x$data %||% list()
    if (length(traces_i)) {
      for (j in seq_along(traces_i)) {
        traces_i[[j]]$xaxis <- "x"
        traces_i[[j]]$yaxis <- ytrace_name
      }
      idx <- trace_index:(trace_index + length(traces_i) - 1L)
      out$x$data[idx] <- traces_i
      trace_index <- trace_index + length(traces_i)
    }

    attrs_i <- p$x$attrs %||% list()
    if (length(attrs_i)) {
      idx <- attrs_index:(attrs_index + length(attrs_i) - 1L)
      out$x$attrs[idx] <- attrs_i
      attrs_index <- attrs_index + length(attrs_i)
    }

    shapes_i <- layout_i$shapes %||% list()
    if (length(shapes_i)) {
      for (j in seq_along(shapes_i)) {
        shape <- shapes_i[[j]]
        shape$xref <- remap_axis_ref(shape$xref %||% "x", "x", "x")
        shape$yref <- remap_axis_ref(shape$yref %||% "y", ytrace_name, "y")
        out$x$layout$shapes[[shape_index]] <- shape
        shape_index <- shape_index + 1L
      }
    }

    annotations_i <- layout_i$annotations %||% list()
    if (length(annotations_i)) {
      for (j in seq_along(annotations_i)) {
        annotation <- annotations_i[[j]]
        annotation$xref <- remap_axis_ref(annotation$xref %||% "x", "x", "x")
        annotation$yref <- remap_axis_ref(annotation$yref %||% "y", ytrace_name, "y")
        out$x$layout$annotations[[annotation_index]] <- annotation
        annotation_index <- annotation_index + 1L
      }
    }

    images_i <- layout_i$images %||% list()
    if (length(images_i)) {
      for (j in seq_along(images_i)) {
        image <- images_i[[j]]
        image$xref <- remap_axis_ref(image$xref %||% "x", "x", "x")
        image$yref <- remap_axis_ref(image$yref %||% "y", ytrace_name, "y")
        out$x$layout$images[[image_index]] <- image
        image_index <- image_index + 1L
      }
    }
  }

  render_hooks <- unlist(lapply(built, function(p) p$jsHooks$render %||%
                                  list()), recursive = FALSE)
  if (length(render_hooks)) out$jsHooks$render <- render_hooks

  deps <- unlist(lapply(built, function(p) p$dependencies %||% list()),
                 recursive = FALSE)
  if (length(deps)) {
    dep_keys <- vapply(
      deps,
      function(d) paste(d$name %||% "", d$version %||% "", d$src$file %||%
                          d$src$href %||% "", sep = "::"),
      character(1)
    )
    out$dependencies <- deps[!duplicated(dep_keys)]
  }

  out
}
