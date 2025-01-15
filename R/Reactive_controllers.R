#' Main controller for browser settings
#' Takes shiny input and converts to a proper list of settings
#' @noRd
click_plot_browser_main_controller <- function(input, tx, cds, libs, df) {
  {
    print(paste("here is gene!", isolate(input$gene)))
    print(paste("here is tx!", isolate(input$tx)))
    display_region <- observed_tx_annotation(isolate(input$tx), tx)
    tx_annotation <- observed_cds_annotation(isolate(input$tx), tx,
                                             isolate(input$other_tx))
    cds_annotation <- observed_cds_annotation(isolate(input$tx), cds,
                                              isolate(input$other_tx))
    uorf_annotation <- observed_uorf_annotation(isolate(input$tx), df,
           isolate(input$other_tx), isolate(input$add_uorfs))
    translon_annotation <- observed_translon_annotation(isolate(input$tx), df(),
                                                isolate(input$other_tx), isolate(input$add_translon))
    customRegions <- c(uorf_annotation, translon_annotation)
    display_region <- genomic_string_to_grl(isolate(input$genomic_region), display_region,
                                            max_size = 1e6, isolate(input$viewMode),
                                            isolate(input$extendLeaders),
                                            isolate(input$extendTrailers))
    dff <- observed_exp_subset(isolate(input$library), libs, df)


    if (isolate(input$withFrames)) {
      withFrames <- libraryTypes(dff, uniqueTypes = FALSE) %in% c("RFP", "RPF", "LSU", "TI")
    } else withFrames <- rep(FALSE, nrow(dff))

    # Hash strings for cache
    hash_strings <- hash_strings_browser(input, dff)

    reads <- try(filepath(dff, "bigwig", suffix_stem = c("_pshifted", "")))
    invalid_reads <- is(reads, "try-error") ||
      (!all(file.exists(unlist(reads, use.names = FALSE))) |
         any(duplicated(unlist(reads, use.names = FALSE))))
    if (invalid_reads) {
      reads <- filepath(dff, "bigwig", suffix_stem = c("_pshifted", ""),
                        base_folders = libFolder(dff, "all"))
    }
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
                   tx_annotation = tx_annotation,
                   reads = reads,
                   custom_sequence = input$customSequence,
                   log_scale = input$log_scale,
                   phyloP = input$phyloP,
                   withFrames = withFrames,
                   hash_bottom = hash_strings[["hash_bottom"]],
                   hash_browser = hash_strings[["hash_browser"]],
                   hash_expression = hash_strings[["hash_expression"]])
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

click_plot_browser_allsamp_controller <- function(input, df, gene_name_list) {
  {
    # browser()
    print(paste("here is gene!", isolate(input$gene)))
    print(paste("here is tx!", isolate(input$tx)))
    id <- isolate(input$tx)
    dff <- df()

    table_path <- collection_path_from_exp(dff, id, isolate(gene_name_list()))
    lib_sizes <- file.path(QCfolder(dff), "totalCounts_mrna.rds")
    if (!file.exists(lib_sizes))
      stop("Count table library size files are not created, missing file totalCounts_mrna.rds",
           " see vignette for more information on how to make these.")
    metadata_field <- isolate(input$metadata)
    clusters <- isolate(input$clusters)
    normalization <- isolate(input$normalization)
    kmer <- isolate(input$kmer)
    min_count <- isolate(input$min_count)
    region_type <- isolate(input$region_type)
    other_gene <- isolate(input$other_gene)
    frame <- isolate(input$frame)
    summary_track <- isolate(input$summary_track)
    ratio_interval <- isolate(input$ratio_interval)
    if (!is.null(ratio_interval)) {
      if (ratio_interval != "") {
        split_ratio <- unlist(strsplit(ratio_interval, ";"))

        temp_interval <- strsplit(split_ratio, ":|-")
        single_input <- lengths(temp_interval) == 1
        if (any(lengths(temp_interval) > 2)) stop("Malformed sort in interval input!")
        if(any(single_input)) {
          temp_interval[which(single_input)] <- lapply(temp_interval[which(single_input)], function(x) rep(x, 2))
        }
        temp_interval <- as.numeric(unlist(temp_interval))

        if (anyNA(temp_interval)) stop("Malformed sort in interval input!")
        ratio_interval <- temp_interval
      } else ratio_interval <- NULL
    }
    subset <- subset_fst_coord_by_region(dff, id, region_type)

    if (!is.null(other_gene) && other_gene != "") {
      print(paste("Sorting on", other_gene))
      other_tx <- tx_from_gene_list(isolate(gene_name_list()), other_gene)[1]
      other_tx_hash <- paste0("sortOther:", other_tx)
    } else {
      other_tx_hash <- other_tx <- NULL
    }

    table_hash <- paste(name(dff), table_path, lib_sizes, clusters, min_count,
                        region_type, metadata_field, normalization, frame,
                        kmer, other_tx_hash, paste(ratio_interval, collapse = ":"),
                        summary_track, sep = "|_|")

    print(paste("Table hash: ", table_hash))
    reactiveValues(dff = dff,
                   id = id,
                   table_path = table_path,
                   lib_sizes = lib_sizes,
                   table_hash = table_hash,
                   metadata_field = metadata_field,
                   normalization = normalization,
                   kmer = kmer,
                   min_count = min_count,
                   subset = subset,
                   region_type = region_type,
                   group_on_tx_tpm = other_tx,
                   ratio_interval = ratio_interval,
                   frame = frame,
                   summary_track = summary_track)
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
    shift_table <- shifts_load(df())
    shift_table <- shift_table[[which(isolate(libs()) == isolate(input$library))]]
    shift_table <- shift_table[fraction %between% c(input$readlength_min, input$readlength_max)]
    if (!is.null(input$p_shifted)) {
      if (!input$p_shifted) additional_extension <- max(abs(shift_table$offsets_start))
    }
  } else {
    warning("Shift table not found!")
    shift_table <- data.table()
  }

  hash_string <- paste(input$extendLeaders + additional_extension,
                       input$extendTrailers + additional_extension,
                       input$normalization, input$region,
                       paste(ORFik:::name_decider(dff, naming = "full"), collapse = "|__|"),
                       input$readlength_min,
                       input$readlength_max, sep = "|__|")

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
                 hash_string = hash_string,
                 reads = reads)
}

click_plot_codon_main_controller <- function(input, tx, cds, libs, df,
                                             length_table) {
  cds_display <- observed_cds_heatmap(isolate(input$tx), cds, length_table,
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
  design <- design(dff, batch.correction.design = TRUE, multi.factor = FALSE)
  target.contrast <- design[1]

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
                 design = design,
                 target.contrast = target.contrast,
                 pairs = pairs, pval = pval,
                 diff_method = diff_method,
                 full = full,
                 hash_string_pre = hash_string_pre,
                 hash_string_full = hash_string_full)
}
