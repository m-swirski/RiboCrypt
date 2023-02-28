click_plot_browser_main_controller <- function(input, tx, cds, libs, df) {
  {
    # browser()
    print(paste("here is gene!", isolate(input$gene)))
    print(paste("here is tx!", isolate(input$tx)))
    display_region <- observed_tx_annotation(isolate(input$tx), tx)
    cds_annotation <- observed_cds_annotation(isolate(input$tx), cds,
                                              isolate(input$other_tx))
    uorf_annotation <- observed_uorf_annotation(names(tx()),
                                                isolate(input$add_uorfs))
    dff <- observed_exp_subset(isolate(input$library), libs, df)
    # customRegions <- load_custom_regions(isolate(input$useCustomRegions), df)
    customRegions <- uorf_annotation
    #reads <- load_reads(dff, "cov")
    reads <- filepath(dff, "bigwig")
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

click_plot_heatmap_main_controller <- function(input, tx, cds, libs, df,
                                               length_table, minFiveUTR = NULL) {
  display_region <- observed_gene_heatmap(isolate(input$tx), tx)
  cds_display <- observed_cds_heatmap(isolate(input$tx), cds, length_table,
                              minFiveUTR = minFiveUTR)
  dff <- observed_exp_subset(isolate(input$library), libs, df)

  additional_extension <- 0
  shift_path <- file.path(ORFik::libFolder(df()), "pshifted", "shifting_table.rds")
  if (file.exists(shift_path)) {
    shift_table <- shifts.load(df())
    shift_table <- shift_table[[which(isolate(libs()) == isolate(input$library))]]
    shift_table <- shift_table[fraction %between% c(input$readlength_min, input$readlength_max)]
    if (!input$p_shifted) additional_extension <- max(abs(shift_table$offsets_start))
  } else {
    warning("Shift table not found!")
    shift_table <- data.table()
  }

  time_before <- Sys.time()
  reads <- load_reads(dff, "covl")

  cat("Library loading: "); print(round(Sys.time() - time_before, 2))
  message("-- Data loading complete")
  reactiveValues(dff = dff,
                 display_region = display_region,
                 extendTrailers = input$extendTrailers + additional_extension,
                 extendLeaders = input$extendLeaders + additional_extension,
                 viewMode = input$viewMode,
                 cds_display = cds_display,
                 region = input$region,
                 readlength_min = input$readlength_min,
                 readlength_max = input$readlength_max,
                 normalization = input$normalization,
                 summary_track = input$summary_track,
                 p_shifted = input$p_shifted,
                 shift_table = shift_table,
                 reads = reads)
}

click_plot_codon_main_controller <- function(input, tx, cds, libs, df,
                                             length_table) {
  cds_display <- observed_cds_heatmap(isolate(input$tx),cds, length_table,
                                      minFiveUTR = 3)
  dff <- observed_exp_subset(isolate(input$library), libs, df)

  time_before <- Sys.time()
  reads <- load_reads(dff, "cov")
  names(reads) <- ORFik:::name_decider(dff, "full")
  cat("Library loading: "); print(round(Sys.time() - time_before, 2))
  message("-- Data loading complete")
  reactiveValues(dff = dff,
                 cds_display = cds_display,
                 reads = reads,
                 codon_score = input$codon_score,
                 filter_value = input$codon_filter_value)
}

click_plot_DEG_main_controller <- function(input, df) {

  dff <- df()[which(df()$condition %in% isolate(input$condition)),]

  time_before <- Sys.time()
  print("experiment subsetting based on condition")
  reactiveValues(dff = dff)
}
