multiOmicsPlot_bottom_panels <- function(reference_sequence, display_range, annotation,
                                         start_codons, stop_codons, custom_motif,
                                         custom_regions, viewMode,
                                         tx_annotation = NULL, collapse_intron_flank = 100) {
  force(display_range)
  # Get sequence and create basic seq panel
  target_seq <- extractTranscriptSeqs(reference_sequence, display_range)
  seq_panel_hits <- createSeqPanelPattern(target_seq[[1]], start_codons = start_codons,
                                          stop_codons = stop_codons, custom_motif = custom_motif)
  seq_panel <- plotSeqPanel(seq_panel_hits, target_seq[[1]])
  # Get the panel for the annotation track
  gene_model_panel <- createGeneModelPanel(display_range, annotation, frame = 1,
                                           tx_annotation = tx_annotation,
                                           custom_regions = custom_regions,
                                           viewMode = viewMode, collapse_intron_flank)
  lines <- gene_model_panel[[2]]
  layers <- max(gene_model_panel[[1]]$layers)
  gene_model_panel <- geneModelPanelPlot(gene_model_panel[[1]])
  return(list(seq_panel = seq_panel, gene_model_panel = gene_model_panel,
              lines = lines, target_seq = target_seq, annotation_layers = layers))
}

multiOmicsPlot_all_track_plots <- function(profiles, withFrames, colors, ylabels,
                                           ylabels_full_name, lines,
                                           frames_type, total_libs,
                                           summary_track, summary_track_type,
                                           BPPARAM) {
  total_libs <- length(profiles)
  
  force(colors)
  force(lines)
  force(ylabels)

  if (frames_type == "animate") {
    plots <- list(getPlotAnimate(rbindlist(profiles, idcol = "file"), withFrames = withFrames[1],
                                colors = colors[1], ylabels = ylabels[1], lines = lines))
  } else {
    if (is(BPPARAM, "SerialParam")) {
      plots <- mapply(function(x,y,z,c,d) createSinglePlot(x,y,z,c,d, lines, type = frames_type, total_libs),
                      profiles, withFrames, colors, ylabels, ylabels_full_name,
                      SIMPLIFY = FALSE)
    } else {
      plots <- bpmapply(function(x,y,z,c,d) createSinglePlot(x,y,z,c,d, lines, type = frames_type, total_libs),
                        profiles, withFrames, colors, ylabels, ylabels_full_name,
                        SIMPLIFY = FALSE, BPPARAM = BPPARAM)
    }
  }

  nplots <- ifelse(frames_type == "animate", 1, length(plots))
  if (summary_track) {
    nplots <- nplots + 1
    plots <- make_summary_track(profiles, plots, withFrames, colors,
                                lines, summary_track_type, nplots)
  }
  return(list(plots = plots, nplots = nplots))
}

multiOmicsPlot_all_profiles <- function(display_range, reads, runs, collection_path, kmers,
                                        kmers_type, frames_type, frames_subset,
                                        withFrames, log_scale, BPPARAM, normalization) {
  force(display_range)
  force(reads)
  force(kmers)
  force(kmers_type)
  force(frames_type)
  force(withFrames)
  force(log_scale)
  force(frames_subset)

  
  if (is(BPPARAM, "SerialParam")) {
    profiles <- mapply(function(x,y,z,b,c,l) getProfileWrapper(display_range, x, y, z, b, c, l, kmers_type,
                                                           type = frames_type, frames_subset = frames_subset, normalization = normalization),
                       reads, runs, collection_path, withFrames, kmers, log_scale, SIMPLIFY = FALSE)
  } else {
    profiles <- bpmapply(function(x,y,z,b,c,l) getProfileWrapper(display_range, x,y,z,b,c,l, kmers_type,
                                                             type = frames_type, frames_subset = frames_subset, normalization = normalization),
                         reads, runs, collection_path, withFrames, kmers, log_scale, SIMPLIFY = FALSE,
                         BPPARAM = BPPARAM)
  }
  return(profiles)
}

multiOmicsPlot_complete_plot <- function(track_panel, bottom_panel, display_range,
                                         proportions, seq_render_dist,
                                         display_sequence, display_dist,
                                         aa_letter_code, input_id, plot_name,
                                         plot_title,  width, height, export.format,
                                         zoom_range = NULL) {
  print("Merging bottom and coverage tracks")
  nplots <- track_panel$nplots
  plots <- browser_plots_highlighted(track_panel$plots, zoom_range)

  gene_model_panel <- bottom_panel$gene_model_panel
  seq_panel <- bottom_panel$seq_panel
  custom_seq_panel <- bottom_panel$custom_bigwig_panels
  without_sequence_track <- display_sequence %in% c("none", FALSE)

  if (without_sequence_track) { # plotly subplot without sequence track
    nplots <- nplots + 2
    plots <- c(plots, list(automateTicksGMP(gene_model_panel), automateTicksX(seq_panel)))
  } else { # plotly subplot with sequence track
    nplots <- nplots + 3
    plots <- c(plots, list(automateTicks(nt_area_template()),
                           automateTicksGMP(gene_model_panel),
                           automateTicksX(seq_panel)))
  }

  if (length(custom_seq_panel) > 0) {
    plots <- c(plots, lapply(custom_seq_panel, automateTicksX))
    nplots_all <- nplots + length(custom_seq_panel)
  } else nplots_all <- nplots

  plots <- lapply(plots, function(x) x  %>% layout(xaxis = list(title = list(font = list(size = 22)), tickfont = list(size = 16)),
                                                   yaxis = list(title = list(font = list(size = 22)), tickfont = list(size = 16),
                                                                rangemode = "tozero", tick0 = 0, dtick = 1,
                                                                zeroline = TRUE)))
  multiomics_plot <- subplot(plots,
                             margin = 0,
                             nrows = nplots_all,
                             heights = proportions,
                             shareX = TRUE,
                             titleY = TRUE, titleX = TRUE)
  if (!without_sequence_track) {
    multiomics_plot <- addJSrender(multiomics_plot, bottom_panel$target_seq,
                                   nplots - 3, seq_render_dist,
                                   display_dist, aa_letter_code, input_id)
  }
  filename <- ifelse(plot_name == "default", names(display_range), plot_name)
  multiomics_plot <- addToImageButtonOptions(multiomics_plot, filename,
                                             width, height, format = export.format)
  if (!is.null(plot_title)) multiomics_plot <- multiomics_plot %>%
    plotly::layout(title = plot_title)
  if (!is.null(zoom_range) && length(zoom_range) == 2) {
    multiomics_plot <- multiomics_plot %>%
      plotly::layout(xaxis = list(range = zoom_range))
  }
  multiomics_plot <- lineDeSimplify(multiomics_plot)

  return(multiomics_plot)
}

#' @importFrom Biostrings nchar translate
#'
multiOmicsPlot_internal <- function(display_range, df, annotation = "cds", reference_sequence = findFa(df),
                                    reads = outputLibs(df, type = "pshifted", output.mode = "envirlist",
                                                       naming = "full", BPPARAM = BiocParallel::SerialParam()),
                                    runs = df[["Run"]],
                                    collection_path = "",
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
                                    summary_track_type = frames_type, export.format = "svg", frames_subset = "all", normalization = normalizations("metabrowser")[1]) {

  multiOmicsController()
  # Get Bottom annotation and sequence panels
  bottom_panel <- multiOmicsPlot_bottom_panels(reference_sequence, display_range, annotation,
                                               start_codons, stop_codons, custom_motif,
                                               custom_regions, viewMode)
  multiOmicsControllerView()

  # Get NGS data track panels
  profiles <- multiOmicsPlot_all_profiles(display_range, reads, runs, collection_path, kmers,
                                          kmers_type, frames_type, frames_subset,
                                          withFrames, log_scale, BPPARAM, normalization)

  track_panel <- multiOmicsPlot_all_track_plots(profiles, withFrames, colors, ylabels,
                                          ylabels_full_name, bottom_panel$lines,
                                          frames_type, total_libs,
                                          summary_track, summary_track_type,
                                          BPPARAM)
  return(multiOmicsPlot_complete_plot(track_panel, bottom_panel, display_range,
                                      proportions, seq_render_dist,
                                      display_sequence, display_dist,
                                      aa_letter_code, input_id, plot_name,
                                      plot_title,  width, height, export.format))
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
                           viewMode, leader_extension, trailer_extension) {
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
          zoom_range <- c(max(as.numeric(start(ir))[1] - 10, 1),
                          min(as.numeric(end(ir))[1] + 10, widthPerGroup(display_range_zoom, FALSE)))
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

hash_strings_browser <- function(
    dff,
    selectedTx, 
    otherTx, 
    addUorfs,
    addTranslons,
    collapsedIntronsWidth,
    genomicRegion,
    extendLeaders,
    extendTrailers,
    zoomRange,
    viewMode,
    withFrames,
    exportFormat,
    summaryTrack,
    summaryTrackType,
    kmer,
    framesType,
    framesSubset,
    customSequence,
    logScale,
    phyloP,
    mapability,
    expressionPlot
    ) {
  full_names <- ORFik:::name_decider(dff, naming = "full")
  hash_bottom <- paste(selectedTx, otherTx,
                       addUorfs,  addTranslons,
                       extendTrailers, extendLeaders,
                       genomicRegion, viewMode,
                       collapsedIntronsWidth,
                       customSequence, phyloP, mapability,
                       collapse = "|_|")
  # Until plot and coverage is split (bottom must be part of browser hash)
  hash_browser <- paste(hash_bottom,
                        full_names,
                        exportFormat,
                        summaryTrack, summaryTrackType,
                        kmer, framesType, withFrames,
                        logScale, zoomRange, framesSubset,
                        collapse = "|_|")
  hash_expression <- paste(full_names, selectedTx,
                           expressionPlot, extendTrailers,
                           extendLeaders, collapse = "|_|")
  hash_strings <- list(hash_bottom = hash_bottom, hash_browser = hash_browser,
                       hash_expression = hash_expression)
  stopifnot(all(lengths(hash_strings) == 1))
  return(hash_strings)
}

