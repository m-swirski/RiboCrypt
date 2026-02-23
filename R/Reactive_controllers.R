#' Main controller for browser settings
#' Takes shiny input and converts to a proper list of settings
#' @noRd
click_plot_browser_main_controller <- function(input, tx, cds, libs, df, gg_theme, user_info) {
  {
    time_before <- controller_init(input, id = "Browser")
    # Annotation
    display_region <- observed_tx_annotation(isolate(input$tx), tx)
    tx_annotation <- observed_cds_annotation(isolate(input$tx), tx,
                                             isolate(input$other_tx))
    cds_annotation <- observed_cds_annotation(isolate(input$tx), cds,
                                              isolate(input$other_tx))
    uorf_annotation <- observed_uorf_annotation(isolate(input$tx), df,
                                                isolate(input$other_tx), isolate(input$add_uorfs))
    translon_annotation <- observed_translon_annotation(isolate(input$tx), df(),
                                                        isolate(input$other_tx), isolate(input$add_translon),
                                                        isolate(input$add_translons_transcode))
    customRegions <- c(uorf_annotation, translon_annotation)

    # View controller
    collapsed_introns_width <- input$collapsed_introns_width
    if (!input$collapsed_introns) collapsed_introns_width <- 0
    if (collapsed_introns_width > 0) {
      tx_annotation <- tx_annotation[tx_annotation %over% flankPerGroup(display_region)]
      display_region_gr <- reduce(unlistGrl(tx_annotation))
      display_region <- groupGRangesBy(display_region_gr, rep(names(display_region), length(display_region_gr)))
    }
    display_region <- genomic_string_to_grl(isolate(input$genomic_region), display_region,
                                            max_size = 1e6, isolate(input$viewMode),
                                            isolate(input$extendLeaders),
                                            isolate(input$extendTrailers),
                                            collapsed_introns_width)
    zoom_range <- get_zoom_range(isolate(input$zoom_range), display_region,
                                 max_size = 1e6, isolate(input$viewMode),
                                 isolate(input$extendLeaders),
                                 isolate(input$extendTrailers))

    if (!is.null(attr(zoom_range, "message"))) {
      showModal(modalDialog(
        title = "Invalid zoom range",
        attr(zoom_range, "message")
      ))
    }

    dff <- observed_exp_subset(isolate(input$library), libs, df)
    if (nrow(dff) > 200) stop("Browser only supports up to 200 libraries for now, use megabrowser!")
    if (isolate(input$withFrames)) {
      withFrames <- libraryTypes(dff, uniqueTypes = FALSE) %in% c("RFP", "RPF", "LSU", "TI")
    } else withFrames <- rep(FALSE, nrow(dff))

    # Hash strings for cache
    hash_strings <- hash_strings_browser(input, dff, collapsed_introns_width)

    if (input$unique_align) uniqueMappers(dff) <- TRUE
    reads <- get_track_paths(dff)
    frames_subset <- input$frames_subset
    use_all_frames <- length(frames_subset) == 0 || any(c("","all") %in% frames_subset)
    if (use_all_frames) frames_subset <- "all"
    frame_colors <- isolate(input$colors)
    colors <- NULL


    cat("-- Browser controller done: "); print(round(Sys.time() - time_before, 2))
    reactiveValues(dff = dff,
                   display_region = display_region,
                   customRegions = customRegions,
                   extendTrailers = input$extendTrailers,
                   extendLeaders = input$extendLeaders,
                   export_format = input$plot_export_format,
                   summary_track = input$summary_track,
                   summary_track_type = input$summary_track_type,
                   viewMode = input$viewMode,
                   collapsed_introns_width = collapsed_introns_width,
                   kmerLength = input$kmer,
                   frames_type = input$frames_type,
                   annotation = cds_annotation,
                   tx_annotation = tx_annotation,
                   reads = reads,
                   custom_sequence = input$customSequence,
                   log_scale = input$log_scale,
                   phyloP = input$phyloP,
                   withFrames = withFrames,
                   zoom_range = zoom_range,
                   frames_subset = frames_subset,
                   mapability = input$mapability,
                   frame_colors = frame_colors,
                   colors = colors,
                   gg_theme = gg_theme,
                   is_cellphone = user_info()$is_cellphone,
                   user_browser_width = user_info()$width,
                   hash_bottom = hash_strings[["hash_bottom"]],
                   hash_browser = hash_strings[["hash_browser"]],
                   hash_expression = hash_strings[["hash_expression"]])
  }
}

click_plot_browser_allsamp_controller <- function(input, df, gene_name_list) {
    time_before <- controller_init(input, id = "Mega Browser")

    # Input copies
    dff <- df()

    id <- isolate(input$tx)
    region_type <- isolate(input$region_type)
    motif <- isolate(input$motif)
    leader_extension <- isolate(input$extendLeaders)
    trailer_extension <- isolate(input$extendTrailers)
    display_annot <- isolate(input$display_annot)
    viewMode <- isolate(input$viewMode)


    # Select motif or gene
    if (isTruthy(motif)) {
      table_path <- meta_motif_files(dff)[motif]
      total_positions <- fst::metadata_fst(table_path)$nrOfRows
      tx_annotation <- display_region <- GRangesList(GRanges(motif, IRanges(1, total_positions), "+"))
      center <- widthPerGroup(tx_annotation, FALSE) / 2
      cds_annotation <- GRangesList(GRanges(motif, IRanges(center - 8, center + 8), "+"))
      names(tx_annotation) <- names(display_region) <- names(cds_annotation) <- motif
      collapsed_introns_width <- 0
      message("Using motif: ", table_path)
    } else {
      # annotation extractors
      annotation_list <- subset_tx_by_region(dff, id, region_type)
      tx_annotation <- display_region <- annotation_list$region
      cds_annotation <- observed_cds_annotation_internal(id,
                                                         annotation_list$cds_annotation,
                                                         isolate(input$other_tx))

      collapsed_introns_width <- input$collapsed_introns_width
      if (!input$collapsed_introns) collapsed_introns_width <- 0
      if (collapsed_introns_width > 0) {
        tx_annotation <- tx_annotation[tx_annotation %over% flankPerGroup(display_region)]
        display_region_gr <- reduce(unlistGrl(tx_annotation))
        display_region <- groupGRangesBy(display_region_gr, rep(names(display_region), length(display_region_gr)))
      }
      display_region <- genomic_string_to_grl(isolate(input$genomic_region), display_region,
                                              max_size = 1e6, viewMode,
                                              leader_extension, trailer_extension,
                                              collapsed_introns_width)
      if (!is.null(leader_extension) && is.numeric(leader_extension) && leader_extension != 0) {
        if (leader_extension > 5e5) stop("Maximum leader extension is 5e5 nt")
        display_region <- extendLeaders(display_region, leader_extension)
      }

      if (!is.null(trailer_extension) && is.numeric(trailer_extension) &&  trailer_extension != 0) {
        if (trailer_extension > 5e5) stop("Maximum leader extension is 5e5 nt")
        display_region <- extendTrailers(display_region, trailer_extension)
      }

      table_path <- collection_path_from_exp(dff, id, isolate(gene_name_list()),
                                             grl_all = display_region)
    }

    # Generic input copies
    clusters <- isolate(input$clusters)
    ratio_interval <- isolate(input$ratio_interval)
    metadata_field <- isolate(input$metadata)
    other_gene <- isolate(input$other_gene)
    enrichment_term <- isolate(input$enrichment_term)
    summary_track <- isolate(input$summary_track)
    normalization <- isolate(input$normalization)
    kmer <- isolate(input$kmer)
    min_count <- isolate(input$min_count)
    frame <- isolate(input$frame)

    # Generic validators and geters
    if (!isTruthy(min_count)) min_count <- 0
    validate_enrichment_term(enrichment_term, clusters, ratio_interval, other_gene,
                             metadata_field)
    ratio_interval <- get_ratio_interval(ratio_interval)
    lib_sizes <- get_lib_sizes_file(dff)
    other_tx_hash <- other_tx <- NULL
    if (isTruthy(other_gene)) {
      print(paste("Sorting on", other_gene))
      other_tx <- tx_from_gene_list(isolate(gene_name_list()), other_gene)[1]
      other_tx_hash <- paste0("sortOther:", other_tx)
    }

    # Hash strings
    table_hash <- paste(name(dff), id, table_path, lib_sizes, clusters, min_count,
                        region_type, paste(metadata_field, collapse = ":"), normalization, frame,
                        kmer, other_tx_hash, paste(ratio_interval, collapse = ":"),
                        leader_extension, trailer_extension,
                        isolate(input$viewMode), collapsed_introns_width, sep = "|_|")
    table_plot_hash <- paste(table_hash, isolate(input$other_tx), input$plotType,
                             summary_track, isolate(input$heatmap_color),
                             isolate(input$color_mult), sep = "|_|")

    timer_done_nice_print("-- Mega Browser controller done: ", time_before)
    reactiveValues(dff = dff,
                   id = id,
                   display_region = display_region,
                   tx_annotation = tx_annotation,
                   annotation = cds_annotation,
                   table_path = table_path,
                   lib_sizes = lib_sizes,
                   metadata_field = metadata_field,
                   normalization = normalization,
                   kmer = kmer,
                   min_count = min_count,
                   viewMode = viewMode,
                   collapsed_introns_width = collapsed_introns_width,
                   subset = NULL,
                   region_type = region_type,
                   group_on_tx_tpm = other_tx,
                   ratio_interval = ratio_interval,
                   frame = frame,
                   clusters = clusters,
                   plotType = isolate(input$plotType),
                   summary_track = summary_track,
                   enrichment_term = enrichment_term,
                   table_hash = table_hash,
                   table_plot_hash = table_plot_hash)
}

controller_init <- function(input, id = "Browser") {
  time_before <- Sys.time()
  print(paste("-",  id, "controller"))
  shinyjs::toggleClass(id = "floating_settings", class = "hidden", condition = TRUE)
  print(paste("Gene:", isolate(input$gene)))
  print(paste("Tx:", isolate(input$tx)))
  if (isolate(input$tx) == "") stop("Empty transcript isoform given!")
  return(time_before)
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

click_plot_heatmap_main_controller <- function(input, tx, cds, libs, df,
                                               minFiveUTR = NULL) {
  display_region <- observed_gene_heatmap(isolate(input$tx), tx)
  cds_display <- observed_cds_heatmap(isolate(input$tx), cds, length_table = NULL,
                              minFiveUTR = minFiveUTR, df = df())
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

  hash_string_anchor <-  paste(input$extendLeaders + additional_extension,
                               input$extendTrailers + additional_extension,
                               input$region, input$customSequence, sep = "|__|")

  hash_string <- paste(hash_string_anchor,
                       input$normalization,
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
                 custom_sequence = input$customSequence,
                 p_shifted = input$p_shifted,
                 shift_table = shift_table,
                 reads = reads,
                 hash_string_anchor = hash_string_anchor,
                 hash_string = hash_string,)
}

click_plot_codon_main_controller <- function(input, tx, cds, libs, df, length_table) {
  cds_display <- observed_cds_heatmap(isolate(input$tx), cds, length_table,
                                      minFiveUTR = 3)
  all_libs <- isolate(input$library)
  background <- if (isTruthy(input$background)) {
    all_libs <- unique(c(all_libs, isolate(input$background)))
    ORFik:::name_decider(observed_exp_subset(isolate(input$background), libs, df), "full")
  }
  dff <- observed_exp_subset(all_libs, libs, df)

  time_before <- Sys.time()
  reads <- load_reads(dff, "cov")
  names <- ORFik:::name_decider(dff, "full")
  names(reads) <- names

  filter_value <- input$codon_filter_value
  normalization <- input$normalization
  differential <- input$differential
  exclude_start_stop <- input$exclude_start_stop
  ratio_thresh <- input$ratio_thresh
  only_significant_difexp <- input$only_significant_difexp


  plot_export_format <- isolate(input$plot_export_format)

  hash_string <- paste(name(dff), names, filter_value, paste(background, collapse = "|lib|"), sep = "|__|")
  hash_string_plot <- paste(hash_string, normalization, differential,
                            exclude_start_stop, ratio_thresh, plot_export_format,
                            only_significant_difexp, sep = "|__|")
  if (differential & length(names) == 1) stop("For differential mode you need at least 2 libraries!")

  cat("Library loading: "); print(round(Sys.time() - time_before, 2))
  message("-- Data loading complete")
  reactiveValues(dff = dff,
                 cds_display = cds_display,
                 reads = reads,
                 normalization = input$normalization,
                 codon_score = input$codon_score,
                 filter_value = filter_value,
                 differential = differential,
                 ratio_thresh = ratio_thresh,
                 exclude_start_stop = exclude_start_stop,
                 background = background,
                 only_significant_difexp = only_significant_difexp,
                 plot_export_format = plot_export_format,
                 hash_string = hash_string,
                 hash_string_plot = hash_string_plot)
}

click_plot_DEG_main_controller <- function(input, df, all_libs, factor = NULL) {
  if (nrow(df()) < 2) stop("Differential expression only allowed for studies with > 1 sample")
  diff_method <- isolate(input$diff_method)
  draw_unregulated <- isolate(input$draw_unnreg)
  conditions <- isolate(input$condition)
  target.contrast <- input$factor
  if (diff_method == "DESeq2" & length(conditions) != 2)
    stop("For DESeq2 method, you must specify 2 levels for the contrast on factor: ", target.contrast)
  # dff <- df()[which(df()$condition %in% conditions),]
  dff <- df()
  design <- factor


  pairs <- list(conditions)
  # pairs <- combn.pairs(unlist(dff[, target.contrast]))
  pval <- isolate(input$pval)

  full <- isolate(input$other_tx)
  libs <- paste(ORFik:::name_decider(dff, naming = "full"), collapse = "")
  group_1 <- input$library1
  group_2 <- input$library2
  plot_export_format <- isolate(input$plot_export_format)

  hash_string_pre <- paste(libs, diff_method, full, paste(group_1, collapse = "|"),
                           paste(group_2, collapse = "|"),
                           sep = "_|-|_")
  hash_string_full <- paste(pval, conditions, plot_export_format,
                            hash_string_pre, sep = "_|-|_")
  hash_string_plot <- paste(draw_unregulated, plot_export_format, hash_string_full,
                            sep = "_|-|_")
  time_before <- Sys.time()
  print("experiment subsetting based on condition")
  reactiveValues(dff = dff, draw_unregulated = draw_unregulated,
                 design = design,
                 target.contrast = target.contrast,
                 all_libs = isolate(all_libs()),
                 group_1 = group_1, group_2 = group_2,
                 pairs = pairs, pval = pval,
                 diff_method = diff_method,
                 full = full,
                 plot_export_format = plot_export_format,
                 hash_string_pre = hash_string_pre,
                 hash_string_full = hash_string_full,
                 hash_string_plot = hash_string_plot)
}
