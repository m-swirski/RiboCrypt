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
  dt[, merged_name := sub(" ",  "-", merged_name, fixed = TRUE)]
  dt[, merged_name := sub("(^-)|(^NA-)",  "", merged_name, perl = TRUE)]
  output_dt <- data.table(value = dt$ensembl_tx_name, label = dt$merged_name)
  if (! is.null(dt$uniprot_id))  output_dt[, uniprot_id := dt$uniprot_id]
  return(output_dt)
}

get_exp <- function(exp_name, experiments, env,
                    exps_dir = ORFik::config()["exp"], page = "") {

  req(exp_name %in% experiments)
  print(paste0("Loading exp: ", exp_name, if (page != "") {paste0("(", page, ")")}))

  exp <- read.experiment(exp_name, output.env = env, validate = FALSE,
                         in.dir = exps_dir)
  print("- New experiment loaded")
  return(exp)
}

bottom_panel_shiny <- function(mainPlotControls) {
  time_before <- Sys.time()
  print("Creating bottom panel..")
  viewMode <- ifelse(mainPlotControls()$viewMode, "genomic", "tx")
  mainPlotControls()$viewMode
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
                                               viewMode,
                                               tx_annotation = mainPlotControls()$tx_annotation,
                                               mainPlotControls()$collapsed_introns_width)
  custom_bigwig_panels <- custom_seq_track_panels(mainPlotControls,
                                                  annotation_list$display_range)
  cat("Done (bottom):"); print(round(Sys.time() - time_before, 2))
  return(c(bottom_panel, annotation_list, custom_bigwig_panels))
}

custom_seq_track_panels <- function(mainPlotControls, display_range) {
  df <- mainPlotControls()$dff
  p <- list()
  if (mainPlotControls()$phyloP) {
    p <- c(p, phylo = list(phylo_custom_seq_track_panel(df, display_range)))
  }
  if (mainPlotControls()$mapability) {
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
  p <- ggplot(data = dt) + geom_bar(aes(y = count, x = position), stat="identity") +
    theme_classic() + theme(axis.title.x = element_blank(),
                            axis.ticks.x = element_blank(),
                            axis.text.x = element_blank(),
                            axis.title.y = element_text(size = rel(0.5)),
                            axis.ticks.y = element_blank(),
                            axis.text.y = element_blank(),
                            plot.margin = unit(c(0,0,0,0), "pt")) +
    scale_x_continuous(expand = c(0,0)) + ylab(ylab)
  return(p)
}


browser_track_panel_shiny <- function(mainPlotControls, bottom_panel, session,
                                      reads = mainPlotControls()$reads,
                                      runs = mainPlotControls()$runs,
                                      collection_path = mainPlotControls()$collection_path,
                                      withFrames = mainPlotControls()$withFrames,
                                      viewMode = ifelse(mainPlotControls()$viewMode, "genomic","tx"),
                                      frames_type = mainPlotControls()$frames_type,
                                      colors = NULL,
                                      kmers = mainPlotControls()$kmerLength,
                                      kmers_type = c("mean", "sum")[1],
                                      ylabels = bamVarName(mainPlotControls()$dff),
                                      lib_to_annotation_proportions = c(0.75,0.25), lib_proportions = NULL,
                                      annotation_proportions = NULL, width = NULL, height = NULL,
                                      plot_name = "default", plot_title = NULL,
                                      display_sequence = "nt", seq_render_dist = 100,
                                      aa_letter_code = c("one_letter", "three_letters")[1],
                                      log_scale = mainPlotControls()$log_scale,
                                      BPPARAM = BiocParallel::SerialParam(),
                                      summary_track = mainPlotControls()$summary_track,
                                      summary_track_type = mainPlotControls()$summary_track_type,
                                      export.format = mainPlotControls()$export_format,
                                      zoom_range = mainPlotControls()$zoom_range,
                                      frames_subset = mainPlotControls()$frames_subset,
                                      normalization = mainPlotControls()$normalization) {
  time_before <- Sys.time()
  print("Creating full browser panel..")
  # Input controller
  multiOmicsControllerView()
  # Get NGS data track panels
  profiles <- multiOmicsPlot_all_profiles(bottom_panel$display_range, reads, runs, collection_path, kmers,
                                          kmers_type, frames_type, frames_subset,
                                          withFrames, log_scale, BPPARAM, normalization)
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
                                       export.format, zoom_range)
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
  return(distribution_plot(boxPlotControls()$dff,
                           boxPlotControls()$display_region,
                           boxPlotControls()$annotation,
                           boxPlotControls()$extendLeaders,
                           boxPlotControls()$extendTrailers))
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
  addResourcePath("tmpuser", dirname(path))  # Ensure the resource path exists
  return(file.path("tmpuser", basename(path)))  # Return only the path
}

click_plot_codon_shiny <- function(mainPlotControls, coverage) {
  click_plot_codon(coverage, min_ratio_change = mainPlotControls$ratio_thresh, min_total_N_codons = 100,
                   exclude_start_stop = mainPlotControls$exclude_start_stop,
                   codon_score = mainPlotControls$codon_score, differential = mainPlotControls$differential)
}


click_plot_codon <- function(dt, min_ratio_change = 1.7, min_total_N_codons = 100,
                             exclude_start_stop = FALSE, codon_score = "percentage",
                             differential = FALSE) {
  message("-- Plotting codon usage")
  if (exclude_start_stop) {
    dt <- dt[grep("(\\*)|(\\#)", seqs, invert = TRUE)]
    dt[, relative_to_max_score := relative_to_max_score / max(relative_to_max_score), by = .(variable, type)]
    dt[, seqs := factor(seqs, levels = levels(seqs)[levels(seqs) %in% unique(seqs)], ordered = TRUE)]
  }

  score_column_name <-
    if (codon_score == "percentage") {
      score_column_name <- "relative_to_max_score"
    } else if (codon_score == "dispersion(NB)") {
      score_column_name <- "dispersion_txNorm"
    } else if (codon_score == "alpha(DMN)") {
      score_column_name <- "alpha"
    } else if (codon_score == "sum") {
      score_column_name <- "sum"
    }

  levels <- which(dt$type == "A" & dt$variable == dt$variable[1])
  levels <- as.character(dt$seqs[levels][order(dt[, ..score_column_name][[1]][levels], decreasing = TRUE)])
  levels <- c(levels, levels(dt$seqs)[!(levels(dt$seqs) %in% levels)])
  dt[, seqs := factor(seqs, levels = levels, ordered = TRUE)]
  dt[, type := factor(type, levels = c("A", "P"), ordered = TRUE)]
  dt <- dt[order(seqs, decreasing = TRUE),][order(type, decreasing = TRUE),]

  score_column <- dt[, ..score_column_name][[1]]


  if (differential) {
    pairs <- ORFik::combn.pairs(unique(dt$variable))
    dt_final <- data.table()
    type <- NULL # avoid BiocCheck error
    for (pair in pairs) {
      sample1 <- dt[variable == pair[1],]
      sample2 <- dt[variable == pair[2],]
      score_column <-
        sample1[, score_column_name, with = FALSE] /
        sample2[, score_column_name, with = FALSE]
      dt_final <- rbindlist(list(dt_final,
                 data.table(variable = paste(sample1$variable,
                                            sample2$variable, sep = " vs "),
                            seqs = sample1$seqs,
                            type = sample1$type,
                            score_column = score_column[[1]],
                            N_sites = paste("S1:", sample1$N.total, "S2:", sample2$N.total),
                            N_min = pmin(sample1$N.total, sample2$N.total))))
    }
    dt <- dt_final
    dt[, simulated_error := sample(c(1,-1))*runif(N_min, 0.1, 0.2)]
    dt[, p_value := ifelse(N_min == 0, NA, wilcox.test(x = rep(score_column + simulated_error, N_min), y = rep(1, N_min), var.equal = FALSE)$p.value),
       by = .(variable, type, seqs)]
    diff_size <- min_ratio_change
    diff_size <- c(diff_size, 1 / diff_size)

    dt[, significant := seq(nrow(dt)) %in% which(N_min >= min_total_N_codons & p_value < 0.05 & (score_column > max(diff_size) | score_column < min(diff_size)))]
    dt[, tooltip_text := paste0(
      "Score: ", round(score_column, 2), "<br>",
      "Sequence: ", seqs, "<br>",
      "Total Sites: ", N_sites, "<br>",
      "P-value:", sprintf("%.2e", p_value)
    )]
    vline_data <- data.frame(x = 1, y = c(min(as.numeric(dt$seqs)), max(as.numeric(dt$seqs))), tooltip_text = "")
    vline_data_p <- data.frame(x = diff_size[1], y = c(min(as.numeric(dt$seqs)), max(as.numeric(dt$seqs))), tooltip_text = "")
    vline_data_m <- data.frame(x = diff_size[2], y = c(min(as.numeric(dt$seqs)), max(as.numeric(dt$seqs))), tooltip_text = "")

    plot <- ggplot(dt, aes(score_column, seqs, text = tooltip_text)) +
      geom_point(color = ifelse(dt$N_min < min_total_N_codons, "red", ifelse(dt$significant == TRUE, "green", "blue"))) +
      geom_line(data = vline_data, aes(x = x, y = y), color = "gray", linetype = "dashed", size = 0.8, alpha = 0.4) +
      geom_line(data = vline_data_p, aes(x = x, y = y), color = "red", linetype = "dashed", size = 0.8, alpha = 0.4) +
      geom_line(data = vline_data_m, aes(x = x, y = y), color = "red", linetype = "dashed", size = 0.8, alpha = 0.4) +
      facet_grid(~ variable + type)
  } else {
    dt[, N_sites := paste(dt$N.total)]
    dt[, N_min := dt$N.tota]
    dt[, tooltip_text := paste0(
      "Score: ", round(score_column, 2), "<br>",
      "Sequence: ", seqs, "<br>",
      "Total Sites: ", N_sites
    )]
    plot <- ggplot(dt, aes(score_column, seqs, text = tooltip_text)) + # fill = score_column
      geom_point(color = ifelse(dt$N_min < 100, "red","blue")) + facet_wrap(variable ~ type, nrow = 1)
    # geom_tile(color = "white") + scale_fill_gradient2(low = "blue", high = "orange", mid = "white")
  }
  plot <- plot + theme_bw() +
    theme(axis.title = element_text(family = "monospace", size = 14, face = "bold"),
          axis.text = element_text(family = "monospace", size = 12),
          strip.background = element_blank(),
          strip.text = element_text(size = 14, face = "bold"),
          strip.placement = "inside",
          plot.margin = margin(t = 20, b = 20, l = 10, r = 10)) +
     ylab("AA:Codon") + xlab("Score (by Ribosome site & library)")

  plotly_plot <- plotly::ggplotly(plot, tooltip = "text")
  return(plotly_plot)
}
