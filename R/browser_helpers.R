multiOmicsPlot_bottom_panels <- function(reference_sequence, display_range, annotation,
                                         start_codons, stop_codons, custom_motif,
                                         custom_regions, viewMode,
                                         tx_annotation = NULL) {
  force(display_range)
  # Get sequence and create basic seq panel
  target_seq <- extractTranscriptSeqs(reference_sequence, display_range)
  seq_panel_hits <- createSeqPanelPattern(target_seq[[1]], start_codons = start_codons,
                                          stop_codons = stop_codons, custom_motif = custom_motif)
  seq_panel <- plotSeqPanel(seq_panel_hits, target_seq[[1]])
  # Get the panel for the annotation track
  gene_model_panel <- createGeneModelPanel(display_range, annotation,
                                           tx_annotation = tx_annotation,
                                           custom_regions = custom_regions,
                                           viewMode = viewMode)
  lines <- gene_model_panel[[2]]
  gene_model_panel <- geneModelPanelPlot(gene_model_panel[[1]])
  return(list(seq_panel = seq_panel, gene_model_panel = gene_model_panel,
              lines = lines, target_seq = target_seq))
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
  if (is(BPPARAM, "SerialParam")) {
    plots <- mapply(function(x,y,z,c,d) createSinglePlot(x,y,z,c,d, lines, type = frames_type, total_libs),
                    profiles, withFrames, colors, ylabels, ylabels_full_name,
                    SIMPLIFY = FALSE)
  } else {
    plots <- bpmapply(function(x,y,z,c,d) createSinglePlot(x,y,z,c,d, lines, type = frames_type, total_libs),
                      profiles, withFrames, colors, ylabels, ylabels_full_name,
                      SIMPLIFY = FALSE, BPPARAM = BPPARAM)
  }
  nplots <- length(plots)
  if (summary_track) {
    nplots <- nplots + 1
    plots <- make_summary_track(profiles, plots, withFrames, colors,
                                lines, summary_track_type, nplots)
  }
  return(list(plots = plots, nplots = nplots))
}

multiOmicsPlot_all_profiles <- function(display_range, reads, kmers,
                                        kmers_type, frames_type,
                                        withFrames, log_scale, BPPARAM) {
  force(display_range)
  force(reads)
  force(kmers)
  force(kmers_type)
  force(frames_type)
  force(withFrames)
  force(log_scale)
  if (is(BPPARAM, "SerialParam")) {
    profiles <- mapply(function(x,y,c,l) getProfileWrapper(display_range, x,y,c,l,kmers_type, type = frames_type),
                       reads, withFrames, kmers, log_scale, SIMPLIFY = FALSE)
  } else {
    profiles <- bpmapply(function(x,y,c,l) getProfileWrapper(display_range, x,y,c,l,kmers_type, type = frames_type),
                         reads, withFrames, kmers, log_scale, SIMPLIFY = FALSE, BPPARAM = BPPARAM)
  }
  return(profiles)
}

multiOmicsPlot_complete_plot <- function(track_panel, bottom_panel, display_range,
                                         proportions, seq_render_dist,
                                         display_sequence, display_dist,
                                         aa_letter_code, input_id, plot_name,
                                         plot_title,  width, height, export.format) {
  nplots <- track_panel$nplots
  plots <- track_panel$plots
  gene_model_panel <- bottom_panel$gene_model_panel
  seq_panel <- bottom_panel$seq_panel
  custom_seq_panel <- bottom_panel$custom_bigwig_panels
  without_sequence_track <- display_sequence %in% c("none", FALSE)
  if (without_sequence_track) { # plotly subplot without sequence track
    nplots <- nplots + 2
    plots <- c(plots, list(automateTicksGMP(gene_model_panel), automateTicksX(seq_panel)))
  } else { # plotly subplot with sequence track
    nplots <- nplots + 3
    plots <- c(plots, list(automateTicks(nt_area_template()), automateTicksGMP(gene_model_panel),
                           automateTicksX(seq_panel)))
  }

  if (!is.null(custom_seq_panel)) {
    plots <- c(plots, list(automateTicksX(custom_seq_panel)))
    nplots_all <- nplots + 1
    proportions <- c(proportions, 0.07)
    proportions <- proportions/sum(proportions)
  } else nplots_all <- nplots

  plots <- lapply(plots, function(x) x  %>% layout(xaxis = list(title = list(font = list(size = 22)), tickfont = list(size = 16)),
                                                   yaxis = list(title = list(font = list(size = 22)), tickfont = list(size = 16))))
  multiomics_plot <- subplot(plots,
                             margin = 0,
                             nrows = nplots_all,
                             heights = proportions,
                             shareX = TRUE,
                             titleY = TRUE, titleX = TRUE)
  if (!without_sequence_track) {
    multiomics_plot <- addJSrender(multiomics_plot, bottom_panel$target_seq,
                                   nplots-3, seq_render_dist,
                                   display_dist, aa_letter_code, input_id)
  }

  filename <- ifelse(plot_name == "default", names(display_range), plot_name)
  multiomics_plot <- addToImageButtonOptions(multiomics_plot, filename,
                                             width, height, format = export.format)
  if (!is.null(plot_title)) multiomics_plot <- multiomics_plot %>% plotly::layout(title = plot_title)
  return(multiomics_plot)
}

#' @importFrom Biostrings nchar translate
#'
multiOmicsPlot_internal <- function(display_range, df, annotation = "cds", reference_sequence = findFa(df),
                                    reads = outputLibs(df, type = "pshifted", output.mode = "envirlist",
                                                       naming = "full", BPPARAM = BiocParallel::SerialParam()),
                                    viewMode = c("tx", "genomic")[1],
                                    custom_regions = NULL,
                                    leader_extension = 0, trailer_extension = 0,
                                    withFrames = libraryTypes(df, uniqueTypes = FALSE) %in% c("RFP", "RPF", "LSU"),
                                    frames_type = "lines", colors = NULL, kmers = NULL, kmers_type = c("mean", "sum")[1],
                                    ylabels = bamVarName(df), lib_to_annotation_proportions = c(0.8,0.2),lib_proportions = NULL,
                                    annotation_proportions = NULL,width = NULL, height = NULL,
                                    plot_name = "default", plot_title = NULL,
                                    display_sequence = c("both","nt", "aa", "none")[1], seq_render_dist = 100,
                                    aa_letter_code = c("one_letter", "three_letters")[1],
                                    annotation_names = NULL, start_codons = "ATG", stop_codons = c("TAA", "TAG", "TGA"),
                                    custom_motif = NULL, log_scale = FALSE, BPPARAM = BiocParallel::SerialParam(),
                                    input_id = "", summary_track = FALSE,
                                    summary_track_type = frames_type, export.format = "svg") {

  multiOmicsController()
  # Get Bottom annotation and sequence panels
  bottom_panel <- multiOmicsPlot_bottom_panels(reference_sequence, display_range, annotation,
                                               start_codons, stop_codons, custom_motif,
                                               custom_regions, viewMode)
  # Get NGS data track panels
  profiles <- multiOmicsPlot_all_profiles(display_range, reads, kmers,
                                          kmers_type, frames_type,
                                          withFrames, log_scale, BPPARAM)
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

genomic_string_to_grl <- function(genomic_string, display_region, max_size = 1e6) {
  input_given <- !is.null(genomic_string) && genomic_string != ""
  if (input_given) {
    gr <- try(GRanges(genomic_string))

    if (is(gr, "GRanges")) {
      if (start(gr) < 1) stop("Position 1 is minimum position to show on a chromosome! (input ", start(gr), ")")
      if (width(gr) > 1e6) stop("Only up to 1 million bases can be shown!")

      try(seqlevelsStyle(gr) <- seqlevelsStyle(display_region)[1], silent = TRUE)

      if (!(as.character(seqnames(gr)) %in% seqnames(seqinfo(df()))))
        stop("Invalid chromosome selected!")
      display_region <- GRangesList(Region = gr)
    } else stop("Malform genomic region: format: chr1:39517672-39523668:+")
  }
  return(display_region)
}

hash_strings_browser <- function(input, dff) {
  full_names <- ORFik:::name_decider(dff, naming = "full")
  hash_bottom <- paste(input$tx, input$other_tx,
                       input$add_uorfs,  input$add_translon,
                       input$extendTrailers, input$extendLeaders,
                       input$genomic_region, input$viewMode,
                       input$customSequence, input$phyloP,
                       collapse = "|_|")
  # Until plot and coverage is split (bottom must be part of browser hash)
  hash_browser <- paste(hash_bottom,
                        full_names,
                        input$plot_export_format,
                        input$summary_track, input$summary_track_type,
                        input$kmer, input$frames_type, input$withFrames,
                        input$log_scale, collapse = "|_|")
  hash_expression <- paste(full_names, input$tx,
                           input$expression_plot, input$extendTrailers,
                           input$extendLeaders, collapse = "|_|")
  hash_strings <- list(hash_bottom = hash_bottom, hash_browser = hash_browser,
                       hash_expression = hash_expression)
  stopifnot(all(lengths(hash_strings) == 1))
  return(hash_strings)
}

