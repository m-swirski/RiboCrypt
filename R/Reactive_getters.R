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
  #print(paste("EXP: ", isolate(dff)))
  return(read.experiment(dff, output.env = env))
}

click_plot_browser <- function(mainPlotControls, session) {
  # browser()
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

click_plot_boxplot <- function(boxPlotControls, session) {
  a <- RiboCrypt:::distribution_plot(boxPlotControls()$dff,
                                     boxPlotControls()$display_region,
                                     boxPlotControls()$annotation,
                                     boxPlotControls()$extendLeaders,
                                     boxPlotControls()$extendTrailers)
  return(a)
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
