get_gene_name_categories <- function(df) {
  print("Updating gene symbol / id table")
  dt <- suppressMessages(symbols(df))
  if (nrow(dt) == 0) {
    # Todo: make it faster
    dt <- mcols(GenomicFeatures::transcripts(loadTxdb(df),
      columns = c("gene_id", "tx_name")
    ))
    dt <- data.table(
      gene_id = unlist(dt$gene_id, use.names = FALSE),
      tx_name = dt$tx_name
    )
    return(data.table(value = dt$tx_name, label = dt$gene_id))
  }
  dt[, merged_name := do.call(paste, .SD, ), .SDcols = c(2, 1)]
  dt[, merged_name := sub(" ", "-", merged_name, fixed = TRUE)]
  dt[, merged_name := sub("(^-)|(^NA-)", "", merged_name, perl = TRUE)]
  output_dt <- data.table(value = dt$ensembl_tx_name, label = dt$merged_name)
  if (!is.null(dt$uniprot_id)) output_dt[, uniprot_id := dt$uniprot_id]
  return(output_dt)
}

get_exp <- function(exp_name, experiments, env,
                    exps_dir = ORFik::config()["exp"], page = "") {
  req(exp_name %in% experiments)
  print(paste0("Loading exp: ", exp_name, if (page != "") {
    paste0("(", page, ")")
  }))

  exp <- read.experiment(exp_name,
    output.env = env, validate = FALSE,
    in.dir = exps_dir
  )
  print("- New experiment loaded")
  return(exp)
}

bottom_panel_shiny <- function(mainPlotControls) {
  time_before <- Sys.time()
  print("Creating bottom panel..")
  viewMode <- ifelse(mainPlotControls()$viewMode, "genomic", "tx")
  df <- mainPlotControls()$dff
  is_cellphone <- mainPlotControls()$is_cellphone

  annotation_list <- annotation_controller(
    df = df,
    display_range = mainPlotControls()$display_region,
    annotation = mainPlotControls()$annotation,
    leader_extension = mainPlotControls()$extendLeaders,
    trailer_extension = mainPlotControls()$extendTrailers,
    viewMode = viewMode
  )

  bottom_panel <- multiOmicsPlot_bottom_panels(
    reference_sequence = findFa(df),
    annotation_list$display_range,
    annotation_list$annotation,
    start_codons = "ATG", stop_codons = c("TAA", "TAG", "TGA"),
    custom_motif = mainPlotControls()$custom_sequence,
    custom_regions = mainPlotControls()$customRegions,
    viewMode,
    tx_annotation = mainPlotControls()$tx_annotation,
    mainPlotControls()$collapsed_introns_width,
    mainPlotControls()$frame_colors,
    mainPlotControls()$gg_theme
  )
  custom_bigwig_panels <- custom_seq_track_panels(
    df, annotation_list$display_range,
    mainPlotControls()$phyloP, mainPlotControls()$mapability
  )
  bottom_panel <- c(bottom_panel, annotation_list, custom_bigwig_panels,
    ncustom = length(custom_bigwig_panels[[1]])
  )

  bottom_panel$bottom_plots <- bottom_plots_to_plotly(bottom_panel, is_cellphone)

  cat("Done (bottom panel):")
  print(round(Sys.time() - time_before, 2))
  return(bottom_panel)
}

bottom_plots_to_plotly <- function(bottom_panel, is_cellphone = FALSE) {
  c(
    list(
      DNA_model = bottom_panel$seq_nt_panel,
      gene_model = automateTicksGMP(bottom_panel$gene_model_panel),
      AA_model = automateTicksAA(bottom_panel$seq_panel, is_cellphone)
    ),
    lapply(bottom_panel$custom_bigwig_panels, automateTicksCustomTrack)
  )
}


custom_seq_track_panels <- function(df, display_range, phyloP, mapability) {
  p <- list()
  if (phyloP) {
    p <- c(p, phylo = list(phylo_custom_seq_track_panel(df, display_range)))
  }
  if (mapability) {
    p2 <- mapability_custom_seq_track_panel(df, display_range)
    p <- c(p, mapability = list(p2))
  }
  return(list(custom_bigwig_panels = p))
}

phylo_custom_seq_track_panel <- function(df, display_range) {
  bw_dir <- file.path(dirname(df@fafile), "phyloP100way")
  if (dir.exists(bw_dir)) {
    bw_track <- list.files(bw_dir, pattern = "\\.phyloP100way\\.bw$", full.names = TRUE)[1]
    if (length(bw_track) == 1) {
      print("- Loading PhyloP track")
      p <- custom_seq_track_panel_bigwig(display_range, bw_track, "P")
      return(p) # If no valid file found
    }
  }
  return(NULL)
}

mapability_custom_seq_track_panel <- function(df, display_range) {
  bw_dir <- file.path(dirname(df@fafile), "mapability")
  if (dir.exists(bw_dir)) {
    bw_track <- list.files(bw_dir, pattern = "28mers_mappability\\.bw$", full.names = TRUE)[1]
    if (length(bw_track) == 1) {
      print("- Loading mapability track")
      p <- custom_seq_track_panel_bigwig(display_range, bw_track, "M")
      return(p) # If no valid file found
    }
  }
  return(NULL)
}

custom_seq_track_panel_bigwig <- function(grl, bigwig_path, ylab) {
  seqlevelsStyle(grl) <- seqlevelsStyle(rtracklayer::BigWigFile(bigwig_path))[1]
  dt <- coveragePerTiling(grl, bigwig_path)
  p <- ggplot(data = dt) +
    geom_bar(aes(y = count, x = position), stat = "identity") +
    theme_classic() +
    theme(
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      axis.title.y = element_text(size = rel(0.5)),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      plot.margin = unit(c(0, 0, 0, 0), "pt")
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    ylab(ylab)
  return(p)
}


browser_track_panel_shiny <- function(mainPlotControls, bottom_panel, session,
                                      reads = mainPlotControls()$reads,
                                      withFrames = mainPlotControls()$withFrames,
                                      viewMode = ifelse(mainPlotControls()$viewMode, "genomic", "tx"),
                                      frames_type = mainPlotControls()$frames_type,
                                      frame_colors = mainPlotControls()$frame_colors,
                                      colors = mainPlotControls()$colors,
                                      kmers = mainPlotControls()$kmerLength,
                                      kmers_type = c("mean", "sum")[1],
                                      ylabels = bamVarName(mainPlotControls()$dff),
                                      lib_to_annotation_proportions = c(0.75, 0.25), lib_proportions = NULL,
                                      annotation_proportions = NULL, width = NULL, height = NULL,
                                      plot_name = "default", plot_title = NULL,
                                      display_sequence = "nt",
                                      seq_render_dist = ifelse(mainPlotControls()$is_cellphone, 100, 250),
                                      aa_letter_code = c("one_letter", "three_letters")[1],
                                      log_scale = mainPlotControls()$log_scale,
                                      BPPARAM = BiocParallel::SerialParam(),
                                      summary_track = mainPlotControls()$summary_track,
                                      summary_track_type = mainPlotControls()$summary_track_type,
                                      export.format = mainPlotControls()$export_format,
                                      zoom_range = mainPlotControls()$zoom_range,
                                      frames_subset = mainPlotControls()$frames_subset) {
  time_before <- Sys.time()
  print("Creating track panel..")
  # Input controller
  multiOmicsControllerView()
  # Get NGS data track panels
  profiles <- multiOmicsPlot_all_profiles(
    bottom_panel$display_range, reads, kmers,
    kmers_type, frames_type, frames_subset,
    withFrames, log_scale, BPPARAM
  )
  track_panel <- multiOmicsPlot_all_track_plots(
    profiles, withFrames, frame_colors,
    colors, ylabels, ylabels_full_name,
    bottom_panel$lines, frames_type,
    summary_track, summary_track_type,
    BPPARAM
  )
  cat("Done (track panel):")
  print(round(Sys.time() - time_before, 2))
  plot <- multiOmicsPlot_complete_plot(track_panel, bottom_panel,
    bottom_panel$display_range,
    proportions, seq_render_dist,
    display_sequence, aa_letter_code,
    input_id = session$ns("selectedRegion"),
    plot_name, plot_title, width, height,
    export.format, zoom_range, frame_colors
  )
  cat("Done (Final panel):")
  print(round(Sys.time() - time_before, 2))
  return(plot)
}

click_plot_browser <- function(mainPlotControls, session) {
  time_before <- Sys.time()
  print("Starting loading + Profile + plot calc")

  a <- RiboCrypt::multiOmicsPlot_ORFikExp(
    display_range = mainPlotControls()$display_region,
    df = mainPlotControls()$dff,
    display_sequence = "nt",
    reads = mainPlotControls()$reads,
    trailer_extension = mainPlotControls()$extendTrailers,
    leader_extension = mainPlotControls()$extendLeaders,
    annotation = mainPlotControls()$annotation,
    viewMode = ifelse(mainPlotControls()$viewMode, "genomic", "tx"),
    kmers = mainPlotControls()$kmerLength,
    frames_type = mainPlotControls()$frames_type,
    custom_regions = mainPlotControls()$customRegions,
    custom_motif = mainPlotControls()$custom_sequence,
    log_scale = mainPlotControls()$log_scale,
    input_id = session$ns("selectedRegion"),
    summary_track = mainPlotControls()$summary_track,
    summary_track_type = mainPlotControls()$summary_track_type,
    export.format = mainPlotControls()$export_format
  )

  cat("lib loading + Coverage calc: ")
  print(round(Sys.time() - time_before, 2))
  return(a)
}

click_plot_boxplot <- function(boxPlotControls, session) {
  return(distribution_plot(
    boxPlotControls()$dff,
    boxPlotControls()$display_region,
    boxPlotControls()$annotation,
    boxPlotControls()$extendLeaders,
    boxPlotControls()$extendTrailers
  ))
}

get_fastq_page <- function(input, libs, df, relative_dir) {
  print("In fastq page")
  dff <- observed_exp_subset(isolate(input$library), libs, df)
  trim_dir <- file.path(libFolder(dff), relative_dir)
  if (!dir.exists(trim_dir)) {
    warning("No valid trim directory")
    return(NULL)
  }
  candidates <- list.files(trim_dir, full.names = TRUE, pattern = "html")
  candidates_base <- gsub("report_", "", sub(".html$", "", basename(candidates)))
  proper_names <- gsub("_Aligned.*", "", ORFik:::remove.file_ext(dff$filepath, basename = TRUE))
  path <- grep(pattern = proper_names, candidates, value = TRUE)

  if (length(path) != 1) {
    hits <- lapply(candidates_base, function(x) grep(x, proper_names))
    path <- candidates[unlist(hits)]
    if (length(path) != 1) {
      warning("No valid html file found in folder!")
      return(NULL)
    }
  }

  print(path)
  addResourcePath("tmpuser", dirname(path)) # Ensure the resource path exists
  return(file.path("tmpuser", basename(path))) # Return only the path
}

click_plot_codon_shiny <- function(mainPlotControls, coverage) {
  click_plot_codon(coverage,
    min_ratio_change = mainPlotControls$ratio_thresh, min_total_N_codons = 100,
    exclude_start_stop = mainPlotControls$exclude_start_stop,
    codon_score = mainPlotControls$codon_score, differential = mainPlotControls$differential,
    background = mainPlotControls$background,
    only_significant_difexp = mainPlotControls$only_significant_difexp,
    format = mainPlotControls$plot_export_format
  )
}

click_plot_codon <- function(dt, min_ratio_change = 1.7, min_total_N_codons = 100,
                             exclude_start_stop = FALSE, codon_score = "percentage",
                             differential = FALSE, background = NULL,
                             only_significant_difexp = FALSE, format = "png") {
  message("-- Plotting codon usage")

  plot <- if (differential) {
    dt <- ORFik::diff_exp_codon(
      dt, min_ratio_change, min_total_N_codons,
      exclude_start_stop, codon_score, background
    )
    ORFik::diff_exp_codon_plot(dt, min_total_N_codons, only_significant_difexp)
  } else {
    ORFik::codon_dotplot(dt, codon_score, min_total_N_codons)
  }

  plotly_plot <- plotly::ggplotly(plot, tooltip = "text") %>%
    plotly::config(toImageButtonOptions = list(format = format, filename = "RC_codon_analysis"), displaylogo = FALSE)
  attr(plotly_plot, "input_data") <- dt
  return(plotly_plot)
}
