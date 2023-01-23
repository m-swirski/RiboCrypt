get_gene_name_categories <- function(df) {
  print("Updating gene symbol / id table")
  dt <- suppressMessages(symbols(df))
  if (nrow(dt) == 0) {
    # Todo: make it faster
    dt <- mcols(GenomicFeatures::transcripts(loadTxdb(df),
                               columns = c("gene_id", "tx_name")))
    dt <- data.table(gene_id = unlist(dt$gene_id, use.names = FALSE),
                     tx_name = dt$tx_name)
    return(data.table(value = dt$tx_name, label = dt$gene_id))
  }
  dt[, merged_name := do.call(paste, .SD, ), .SDcols = c(2,1)]
  dt[, merged_name := gsub(" ",  "-", merged_name)]
  dt[, merged_name := gsub("^-",  "", merged_name)]
  return(data.table(value = dt$ensembl_tx_name, label = dt$merged_name))
}

get_exp <- function(dff, experiments, env) {
  print("testing exp")
  req(dff %in% experiments)
  print("New experiment loaded")
  return(read.experiment(dff, output.env = env))
}

click_plot_browser_main_controller <- function(input, tx, cds, libs, df) {
  {
    print(paste("here is gene!", isolate(input$gene)))
    print(paste("here is tx!", isolate(input$tx)))
    display_region <- observed_tx_annotation(isolate(input$tx), tx)
    cds_annotation <- observed_cds_annotation(isolate(input$tx), cds,
                                              isolate(input$other_tx))
    dff <- observed_exp_subset(isolate(input$library), libs, df)
    customRegions <- load_custom_regions(isolate(input$useCustomRegions), df)

    reads <- load_reads(dff, "cov")
    reactiveValues(dff = dff,
                   display_region = display_region,
                   customRegions = customRegions,
                   extendTrailers = input$extendTrailers,
                   extendLeaders = input$extendLeaders,
                   export_format = input$plot_export_format,
                   summary_track = input$summary_track,
                   summary_track_type = input$summary_track_type,
                   viewMode = input$viewMode,
                   kmerLength = input$kmer,
                   frames_type = input$frames_type,
                   annotation = cds_annotation,
                   reads = reads)
  }
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
    viewMode = ifelse(mainPlotControls()$viewMode, "genomic","tx"),
    kmers = mainPlotControls()$kmerLength,
    frames_type = mainPlotControls()$frames_type,
    custom_regions = mainPlotControls()$customRegions,
    input_id = session$ns("selectedRegion"),
    summary_track = mainPlotControls()$summary_track,
    summary_track_type = mainPlotControls()$summary_track_type,
    export.format = mainPlotControls()$export_format
  )
  cat("lib loading + Coverage calc: "); print(round(Sys.time() - time_before, 2))
  return(a)
}

click_plot_heatmap_main_controller <- function(input, tx, cds, libs, df) {
  display_region <- observed_gene_heatmap(isolate(input$tx), tx)
  cds_display <- observed_cds_heatmap(isolate(input$tx),cds)
  dff <- observed_exp_subset(isolate(input$library), libs, df)


  time_before <- Sys.time()
  reads <- load_reads(dff, "covl")
  cat("Library loading: "); print(round(Sys.time() - time_before, 2))
  message("-- Data loading complete")
  reactiveValues(dff = dff,
                 display_region = display_region,
                 extendTrailers = input$extendTrailers,
                 extendLeaders = input$extendLeaders,
                 viewMode = input$viewMode,
                 cds_display = cds_display,
                 region = input$region,
                 readlength_min = input$readlength_min,
                 readlength_max = input$readlength_max,
                 normalization = input$normalization,
                 reads = reads)
}

click_plot_heatmap_main_controller <- function(input, tx, cds, libs, df) {
  display_region <- observed_gene_heatmap(isolate(input$tx), tx)
  cds_display <- observed_cds_heatmap(isolate(input$tx),cds)
  dff <- observed_exp_subset(isolate(input$library), libs, df)


  time_before <- Sys.time()
  reads <- load_reads(dff, "covl")
  cat("Library loading: "); print(round(Sys.time() - time_before, 2))
  message("-- Data loading complete")
  reactiveValues(dff = dff,
                 display_region = display_region,
                 extendTrailers = input$extendTrailers,
                 extendLeaders = input$extendLeaders,
                 viewMode = input$viewMode,
                 cds_display = cds_display,
                 region = input$region,
                 readlength_min = input$readlength_min,
                 readlength_max = input$readlength_max,
                 normalization = input$normalization,
                 reads = reads)
}

click_plot_codon_main_controller <- function(input, tx, cds, libs, df) {
  cds_display <- observed_cds_heatmap(isolate(input$tx),cds)
  dff <- observed_exp_subset(isolate(input$library), libs, df)

  time_before <- Sys.time()
  reads <- load_reads(dff, "covl")
  names(reads) <- ORFik:::name_decider(dff, "full")
  cat("Library loading: "); print(round(Sys.time() - time_before, 2))
  message("-- Data loading complete")
  reactiveValues(dff = dff,
                 cds_display = cds_display,
                 reads = reads,
                 filter_value = input$codon_filter_value)
}