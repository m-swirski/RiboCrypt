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
  dt[, merged_name := gsub("(^-)|(^NA-)",  "", merged_name)]
  output_dt <- data.table(value = dt$ensembl_tx_name, label = dt$merged_name)
  if (! is.null(dt$uniprot_id))  output_dt$uniprot_id <- dt$uniprot_id
  return(output_dt)
}

get_gene_name_categories_collection <- function(df) {
  valid <- list.files(file.path(resFolder(df), "collection_tables/"))
  valid <- gsub("\\.fst", "", valid)

  all_genes <- get_gene_name_categories(df)
  return(all_genes[value%in% valid,])
}

get_exp <- function(dff, experiments, env) {
  print("testing exp")
  req(dff %in% experiments)
  print("New experiment loaded")
  #print(paste("EXP: ", isolate(dff)))
  return(read.experiment(dff, output.env = env, validate = FALSE))
}

bottom_panel_shiny <- function(mainPlotControls) {
  time_before <- Sys.time()
  print("Creating bottom panel..")
  viewMode <- ifelse(mainPlotControls()$viewMode, "genomic", "tx")
  df <- mainPlotControls()$dff
  annotation_list <- annotation_controller(df = df,
                                           display_range = mainPlotControls()$display_region,
                                           annotation = mainPlotControls()$annotation,
                                           leader_extension = mainPlotControls()$extendLeaders,
                                           trailer_extension = mainPlotControls()$extendTrailers,
                                           viewMode = viewMode)

  bottom_panel <- multiOmicsPlot_bottom_panels(reference_sequence = findFa(df),
                                               annotation_list$display_range,
                                               annotation_list$annotation,
                                               start_codons = "ATG", stop_codons = c("TAA", "TAG", "TGA"),
                                               custom_motif = mainPlotControls()$custom_sequence,
                                               custom_regions = mainPlotControls()$customRegions,
                                               viewMode)
  custom_bigwig_panels <- custom_seq_track_panels(mainPlotControls)
  cat("Done (bottom):"); print(round(Sys.time() - time_before, 2))
  return(c(bottom_panel, annotation_list, custom_bigwig_panels))
}

custom_seq_track_panels <- function(mainPlotControls) {
  if (!mainPlotControls()$phyloP) return(NULL)
  df <- mainPlotControls()$dff
  phylo_dir <- file.path(dirname(df@fafile), "phyloP100way")
  if (dir.exists(phylo_dir)) {
    phylo_track <- list.files(phylo_dir, pattern = "\\.phyloP100way\\.bw$", full.names = TRUE)[1]
    if (length(phylo_track) == 1) {
      print("- Loading PhyloP track")
      grl <- mainPlotControls()$display_region
      if (mainPlotControls()$viewMode) {
        rl <- unlistGrl(flankPerGroup(grl))
      }
      seqlevelsStyle(grl) <- seqlevelsStyle(rtracklayer::BigWigFile(phylo_track))

      rl <- ranges(grl)
      names(rl) <- seqnamesPerGroup(grl, FALSE)

      res <- rtracklayer::import.bw(phylo_track, as = "NumericList",
                                    which = rl)
      dt <- data.table(phyloP = unlist(res, use.names = FALSE))
      dt[, position := seq.int(.N)]
      p <- ggplot(data = dt) + geom_bar(aes(y = phyloP, x = position), stat="identity") +
        theme_classic() + theme(axis.title.x = element_blank(),
                                axis.ticks.x = element_blank(),
                                axis.text.x = element_blank(),
                                axis.title.y = element_blank(),
                                axis.ticks.y = element_blank(),
                                axis.text.y = element_blank(),
                                plot.margin = unit(c(0,0,0,0), "pt")) +
        scale_x_continuous(expand = c(0,0))
      return(list(custom_bigwig_panels = p))
    }
  }
  return(NULL) # If no valid file found
}

browser_track_panel_shiny <- function(mainPlotControls, bottom_panel, session,
                                      reads = mainPlotControls()$reads,
                                      withFrames = mainPlotControls()$withFrames,
                                      viewMode = ifelse(mainPlotControls()$viewMode, "genomic","tx"),
                                      frames_type = mainPlotControls()$frames_type,
                                      colors = NULL,
                                      kmers = mainPlotControls()$kmerLength,
                                      kmers_type = c("mean", "sum")[1],
                                      ylabels = bamVarName(mainPlotControls()$dff),
                                      lib_to_annotation_proportions = c(0.8,0.2), lib_proportions = NULL,
                                      annotation_proportions = NULL, width = NULL, height = NULL,
                                      plot_name = "default", plot_title = NULL,
                                      display_sequence = "nt", seq_render_dist = 100,
                                      aa_letter_code = c("one_letter", "three_letters")[1],
                                      log_scale = mainPlotControls()$log_scale,
                                      BPPARAM = BiocParallel::SerialParam(),
                                      summary_track = mainPlotControls()$summary_track,
                                      summary_track_type = mainPlotControls()$summary_track_type,
                                      export.format = mainPlotControls()$export_format) {
  time_before <- Sys.time()
  print("Creating full browser panel..")
  # Input controller
  multiOmicsControllerView()
  # Get NGS data track panels
  # browser()
  profiles <- multiOmicsPlot_all_profiles(bottom_panel$display_range, reads, kmers,
                                          kmers_type, frames_type,
                                          withFrames, log_scale, BPPARAM)
  track_panel <- multiOmicsPlot_all_track_plots(profiles, withFrames, colors, ylabels,
                                                ylabels_full_name, bottom_panel$lines,
                                                frames_type, total_libs,
                                                summary_track, summary_track_type,
                                                BPPARAM)
  plot <- multiOmicsPlot_complete_plot(track_panel, bottom_panel,
                                       bottom_panel$display_range,
                                       proportions, seq_render_dist,
                                       display_sequence, display_dist,
                                       aa_letter_code,
                                       input_id = session$ns("selectedRegion"),
                                       plot_name, plot_title, width, height,
                                       export.format)
  cat("Done (Full):"); print(round(Sys.time() - time_before, 2))
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
    viewMode = ifelse(mainPlotControls()$viewMode, "genomic","tx"),
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

  cat("lib loading + Coverage calc: "); print(round(Sys.time() - time_before, 2))
  return(a)
}

click_plot_boxplot <- function(boxPlotControls, session) {
  a <- RiboCrypt:::distribution_plot(boxPlotControls()$dff,
                                     boxPlotControls()$display_region,
                                     boxPlotControls()$annotation,
                                     boxPlotControls()$extendLeaders,
                                     boxPlotControls()$extendTrailers)
  return(a)
}

compute_collection_table_shiny <- function(mainPlotControls,
                                      path = mainPlotControls()$table_path,
                                      lib_sizes = mainPlotControls()$lib_sizes,
                                      df = mainPlotControls()$dff,
                                      metadata_field = mainPlotControls()$metadata_field,
                                      normalization = mainPlotControls()$normalization,
                                      kmer = mainPlotControls()$kmer,
                                      min_count = mainPlotControls()$min_count,
                                      subset = mainPlotControls()$subset,
                                      group_on_tx_tpm = mainPlotControls()$group_on_tx_tpm,
                                      metadata) {
  if (is.null(metadata)) stop("Metadata not defined, no metabrowser allowed for now!")
  time_before <- Sys.time()
  cat("Starting loading + Profile + plot calc\n")
  dtable <- compute_collection_table(path, lib_sizes, df, metadata_field,
                                     normalization, kmer, metadata, min_count,
                                     as_list = TRUE, subset = subset,
                                     group_on_tx_tpm = group_on_tx_tpm)
  cat("Done: lib loading + Coverage calc: "); print(round(Sys.time() - time_before, 2))
  return(dtable)
}

allsamples_metadata_clustering <- function(values, plot) {
  time_before <- Sys.time()
  print("Starting metabrowser clustering info")

  pdf(NULL) # TODO: Make a better fix for blank pdf write
  row_orders <- suppressWarnings(ComplexHeatmap::row_order(plot))
  dev.off()
  orders <- unlist(row_orders, use.names = FALSE)
  clustering_was_done <- is.list(row_orders)

  if (!clustering_was_done) {
    orders <- order(values) # Order by variable instead of cluster
  }

  meta <- data.table(grouping = values, order = orders)
  meta <- meta[meta$order, ]
  if (clustering_was_done) {
    meta[, cluster := rep(seq(length(row_orders)), lengths(row_orders))]
  }
  meta[, index := .I]
  cat("metabrowser clustering info done"); print(round(Sys.time() - time_before, 2))
  return(meta)
}

allsamples_sidebar <- function(meta) {
  numeric_grouping <- is.numeric(meta$grouping)

  gg <- ggplot(meta, aes(y = rev(index), x = factor(1), fill = grouping)) +
    theme_void() +
    labs(x = NULL, y = NULL, title = NULL) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(panel.background = element_rect(fill = "transparent", colour = NA),
          plot.background = element_rect(fill = "transparent", colour = NA),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          plot.margin = unit(c(0, 0, 0, 0), "cm"),
          legend.position = "none")

  if (numeric_grouping) {
    meta[, grouping_numeric := grouping]
    meta[, grouping := cut(grouping, breaks=5)]
    gg_tpm <- gg + geom_line(aes(y = grouping_numeric, x = rev(index), fill = NULL)) +
      coord_flip()
    plotly_tpm <- ggplotly(gg_tpm, tooltip="text")
  }

  gg_groups <- gg + geom_raster()
  res <- ggplotly(gg_groups, tooltip="text")

  if (numeric_grouping) {
    res <- subplot(plotly_tpm, res, nrows = 1)
  }
  return(res %>% plotly::config(displayModeBar = FALSE))
}


allsamples_meta_stats_shiny <- function(meta) {
  if (is.numeric(meta$grouping)) {
    meta[, grouping := cut(grouping, breaks=5)]
  }
  dt <- allsamples_meta_stats(meta)
  # Add Chi squared significane coloring
  datatable(round(dt, 2)) %>% formatStyle(columns = seq(ncol(dt)),
     backgroundColor = styleInterval(c(-3, 3), c('yellow', 'white', 'yellow')))
}

allsamples_meta_stats <- function(meta) {
  time_before <- Sys.time()
  print("Starting metabrowser statistics")
  res <- copy(meta)
  res[, index := NULL]
  res[, grouping := as.character(grouping)]
  clustering_was_done <- !is.null(res$cluster)
  if (clustering_was_done) {
    concat_table <- table(res$grouping, res$cluster)
    chi_test <- chisq.test(concat_table)
    res <- chi_test$stdres
  } else {
    concat_table <- table(res$grouping)
    res <- matrix(concat_table, ncol = 1)
    rownames(res) <- names(concat_table)
    colnames(res) <- "counts"
  }
  cat("metabrowser statistics done"); print(round(Sys.time() - time_before, 2))
  return(as.data.frame.matrix(res))
}

get_meta_browser_plot <- function(table, color_theme, clusters = 1,
                                  color_mult = 3) {
  colors <- if (color_theme == "default (White-Blue)") {
    c("white", "lightblue", rep("blue", 4 + color_mult), "navy", "black")
  } else if (color_theme == "Matrix (black,green,red)") {
    c("#000000", "#2CFA1F", "yellow2", rep("#FF2400", color_mult))
  } else stop("Invalid color theme!")
  cat("Creating metabrowser heatmap\n")
  ComplexHeatmap::Heatmap(t(table), show_row_dend = FALSE,
                          cluster_columns = FALSE,
                          cluster_rows = FALSE,
                          use_raster = TRUE,  raster_quality = 5,
                          km = clusters,
                          col =  colors, show_row_names = FALSE)
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
  candidates_base <- gsub("report_", "", gsub(".html$", "", basename(candidates)))
  proper_names <- gsub("_Aligned.*", "", ORFik:::remove.file_ext(dff$filepath,basename = TRUE))
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
  addResourcePath("tmpuser", dirname(path))
  path <- file.path("tmpuser", basename(path))
  page <- tags$iframe(seamless="seamless", src= path, width=1000, height=900)
}

click_plot_codon <- function(input, coverage) {
  message("-- Plotting codon usage")

  score_column <-
    if (input$codon_score == "percentage") {
      score_column_name <- "relative_to_max_score"
      coverage()$relative_to_max_score
    } else if (input$codon_score == "dispersion(NB)") {
      score_column_name <- "dispersion_txNorm"
      coverage()$dispersion_txNorm
    } else if (input$codon_score == "alpha(DMN)") {
      score_column_name <- "alpha"
      coverage()$alpha
    } else if (input$codon_score == "sum") {
      score_column_name <- "sum"
      coverage()$sum
    }

  if (input$differential) {
    pairs <- ORFik::combn.pairs(unique(coverage()$variable))
    dt <- data.table()
    type <- NULL # avoid BiocCheck error
    for (pair in pairs) {
      sample1 <- coverage()[variable == pair[1],]
      sample2 <- coverage()[variable == pair[2],]
      score_column <-
        sample1[, score_column_name, with = FALSE] /
        sample2[, score_column_name, with = FALSE]
      dt <- rbindlist(list(dt,
                 data.table(variable = paste(sample1$variable,
                                            sample2$variable, sep = " vs "),
                            seqs = sample1$seqs,
                            type = sample1$type,
                            score_column = score_column[[1]])))
    }
    plotly::ggplotly(ggplot(dt,
                            aes(score_column, seqs)) +
                       geom_point(color = "blue") +
                       scale_fill_gradient2(low = "blue", high = "orange",
                                            mid = "white") +
                       theme(axis.text.y = element_text(family = "monospace")) +
                       facet_grid(type ~ variable))
  } else {
    plotly::ggplotly(ggplot(coverage(),
                            aes(type, seqs, fill = score_column)) +
                       geom_tile(color = "white") +
                       scale_fill_gradient2(low = "blue", high = "orange",
                                            mid = "white") +
                       theme(axis.text.y = element_text(family = "monospace")) +
                       facet_wrap(coverage()$variable, ncol = 4))
  }

}
