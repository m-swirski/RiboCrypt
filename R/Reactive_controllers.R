click_plot_browser_main_controller <- function(input, tx, cds, libs, df) {
  {
    # browser()
    print(paste("here is gene!", isolate(input$gene)))
    print(paste("here is tx!", isolate(input$tx)))
    display_region <- observed_tx_annotation(isolate(input$tx), tx)
    cds_annotation <- observed_cds_annotation(isolate(input$tx), cds,
                                              isolate(input$other_tx))
    uorf_annotation <- observed_uorf_annotation(isolate(input$tx), df,
           isolate(input$other_tx), isolate(input$add_uorfs))
    if (!is.null(isolate(input$genomic_region)) & isolate(input$genomic_region) != "") {
      gr <- try(GRanges(isolate(input$genomic_region)))

      if (is(gr, "GRanges")) {
        if (width(gr) > 1e6) stop("Only up to 1 million bases can be shown!")
        seqlevelsStyle(gr) <- seqlevelsStyle(display_region)[1]
        if (!(as.character(seqnames(gr)) %in% seqnames(seqinfo(df()))))
          stop("Invalid chromosome selected!")
        display_region <- GRangesList(Region = gr)
      } else stop("Malform genomic region: format: chr1:39517672-39523668:+")
    }
    dff <- observed_exp_subset(isolate(input$library), libs, df)
    # customRegions <- load_custom_regions(isolate(input$useCustomRegions), df)
    customRegions <- uorf_annotation
    #reads <- load_reads(dff, "cov")
    reads <- filepath(dff, "bigwig", suffix_stem = c("_pshifted", ""))
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
                   reads = reads,
                   custom_sequence = input$customSequence)
  }
}

click_plot_browser_new_controller <- function(input, tx, cds, libs, df) {
  {
    print(paste("here is gene!", isolate(input$gene)))
    print(paste("here is tx!", isolate(input$tx)))
    display_region <- observed_tx_annotation(isolate(input$tx), tx)
    annotation <- observed_cds_annotation(isolate(input$tx), cds,
                                              isolate(input$other_tx))
    dff <- observed_exp_subset(isolate(input$library), libs, df)
    customRegions <- load_custom_regions(isolate(input$useCustomRegions), df)

    #reads <- load_reads(dff, "cov")
    reads <- filepath(dff, "bigwig", suffix_stem = c("_pshifted", ""))
    trailer_extension <- input$extendTrailers
    leader_extension <- input$extendLeaders
    export.format <- input$plot_export_format
    summary_track <- input$summary_track
    summary_track_type <- input$summary_track_type
    viewMode <- input$viewMode
    kmers <- input$kmer
    frames_type <- input$frames_type

    #other defaults
    annotation_names <- NULL
    display_sequence <- "both"
    withFrames <- libraryTypes(dff, uniqueTypes = FALSE) %in% c("RFP", "RPF", "LSU")
    lib_proportions <- NULL
    colors = NULL
    ylabels = bamVarName(dff)
    lib_to_annotation_proportions <- c(0.8,0.2)

    multiOmicsController()

    reactiveValues(dff = dff,
                   display_region = display_region,
                   customRegions = customRegions,
                   extend_trailers = extend_trailers,
                   extend_leaders = extend_trailers,
                   export.format = export.format,
                   summary_track = summary_track,
                   summary_track_type = summary_track_type,
                   viewMode = viewMode,
                   kmers = kmers,
                   frames_type = frames_type,
                   annotation = annotation,
                   reads = reads,
                   withFrames = withFrames,
                   )
  }
}
# click_plot_boxplot_main_controller <- function(input, tx, libs, df) {
#   dff <- observed_exp_subset(isolate(input$library), libs, df)
#   customRegions <- load_custom_regions(isolate(input$useCustomRegions), df)
#
#   reactiveValues(dff = dff,
#                  tx  = isolate(input$tx),
#                  customRegions = customRegions,
#                  export_format = input$plot_export_format)
# }

click_plot_browser_allsamp_controller <- function(input, df) {
  {
    # browser()
    print(paste("here is gene!", isolate(input$gene)))
    print(paste("here is tx!", isolate(input$tx)))
    dff <- df()

    table_path <- file.path(resFolder(dff), "collection_tables", paste0(isolate(input$tx), ".fst"))
    if (!file.exists(dirname(table_path)))
      stop("There is no collection fst tables directory for this organism,",
           " see vignette for more information on how to make these.")
    if (!file.exists(table_path)) stop("Gene has no precomputed table, try another one!")
    lib_sizes <- file.path(QCfolder(dff), "totalCounts_mrna.rds")
    if (!file.exists(lib_sizes))
      stop("Count table library size files are not created, missing file totalCounts_mrna.rds",
           " see vignette for more information on how to make these.")
    metadata_field <- isolate(input$metadata)
    clusters <- isolate(input$clusters)
    normalization <- isolate(input$normalization)
    kmer <- isolate(input$kmer)
    table_hash <- paste(name(dff), table_path, lib_sizes, clusters,
                        metadata_field, normalization, kmer, sep = "|_|")

    print(paste("Table hash: ", table_hash))
    reactiveValues(dff = dff,
                   table_path = table_path,
                   lib_sizes = lib_sizes,
                   table_hash = table_hash,
                   metadata_field = metadata_field,
                   normalization = normalization,
                   kmer = kmer)
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
    if (!is.null(input$p_shifted)) {
      if (!input$p_shifted) additional_extension <- max(abs(shift_table$offsets_start))
    }
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
  if (nrow(df()) < 2) stop("Differential expression only allowed for studies with > 1 sample")
  draw_unregulated <- isolate(input$draw_unnreg)
  conditions <- isolate(input$condition)
  # dff <- df()[which(df()$condition %in% conditions),]
  dff <- df()
  target.contrast <- "condition"
  pairs <- combn.pairs(unlist(dff[, target.contrast]))
  pval <- isolate(input$pval)
  diff_method = isolate(input$diff_method)
  full <- isolate(input$other_tx)
  libs <- paste(ORFik:::name_decider(dff, naming = "full"), collapse = "")
  hash_string_pre <- paste(draw_unregulated,
                           libs, diff_method, full,
                           sep = "_|-|_")
  hash_string_full <- paste(pval, conditions,
                            hash_string_pre, sep = "_|-|_")
  time_before <- Sys.time()
  print("experiment subsetting based on condition")
  reactiveValues(dff = dff, draw_unregulated = draw_unregulated,
                 target.contrast = target.contrast,
                 pairs = pairs, pval = pval,
                 diff_method = diff_method,
                 full = full,
                 hash_string_pre = hash_string_pre,
                 hash_string_full = hash_string_full)
}
