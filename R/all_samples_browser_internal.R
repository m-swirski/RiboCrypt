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
                                           split_by_frame = mainPlotControls()$frame,
                                           ratio_interval = mainPlotControls()$ratio_interval,
                                           metadata) {
  if (is.null(metadata)) stop("Metadata not defined, no metabrowser allowed for now!")
  time_before <- Sys.time()
  cat("Starting loading + Profile + plot calc\n")
  dtable <- compute_collection_table(path, lib_sizes, df, metadata_field,
                                     normalization, kmer, metadata, min_count,
                                     as_list = TRUE, subset = subset,
                                     group_on_tx_tpm = group_on_tx_tpm,
                                     split_by_frame = split_by_frame,
                                     ratio_interval = ratio_interval)
  cat("Done: lib loading + Coverage calc: "); print(round(Sys.time() - time_before, 2))
  return(dtable)
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
                          col =  colors, show_row_names = FALSE,
                          show_heatmap_legend = FALSE)
}

#' Full plot for allsamples browser
#'
#' @param m data.table of coverage per sample (wide format)
#' @param heatmap ComplexHeatmap object of plot from 'm'
#' @param id id of transcript
#' @param df ORFik experiment
#' @param sumary logical, default TRUE (add top plot)
#' @param annotation logical, default TRUE (add bottom annotation track)
#' @param region_type character, "what is the coverage region?" Usually full mrna:
#' "mrna" or "leader+cds".
#' @param rel_heights numeric < 1, default: c(0.2, 0.75, 0.05).
#' Relative heights, sum to 1 and must be length 3.
#' @importFrom cowplot plot_grid
get_meta_browser_plot_full <- function(m, heatmap, id, df,
                                       summary = TRUE, annotation = TRUE,
                                       region_type,
                                       rel_heights = c(0.2, 0.75, 0.05)) {
  stopifnot(is.numeric(rel_heights) && length(rel_heights) == 3)
  gene_model_panel <- summary_plot <- NULL
  if (!summary & !annotation) return(heatmap)

  if (summary) {
    summary_plot <- summary_track_allsamples(m)
  }

  if (annotation) {
    gene_model_panel <- annotation_track_allsamples(df, id, region_type,
                                                    tx_width = ncol(heatmap))
  }

  grob <- grid::grid.grabExpr(draw(heatmap))
  to_use_logicals <- c(summary, TRUE, annotation)
  final_plot <- cowplot::plot_grid(plotlist = list(summary_plot, grob, gene_model_panel)[to_use_logicals],
                                   ncol = 1, rel_heights = rel_heights[to_use_logicals])
  return(final_plot)
}

summary_track_allsamples <- function(m, summary_track_type = "area", as_plotly = FALSE) {

  summary_profile <- data.table(count = rowSums(m))
  summary_profile[, `:=`(position = seq.int(.N)) ]
  summary_profile[, `:=`(frame = factor((position-1) %% 3)) ]
  summary_plot <- createSinglePlot(summary_profile, TRUE, 1, "",
                                   FALSE, lines = NULL,
                                   type = summary_track_type,
                                   flip_ylabel = FALSE, as_plotly = as_plotly)
  if (!as_plotly) {
    summary_plot +
      theme(axis.text.y=element_blank(),  axis.ticks.y=element_blank())
  }
  return(summary_plot)
}

annotation_track_allsamples <- function(df, id, region_type, tx_width) {
  lines <- NULL
  withFrames <- TRUE
  colors <- TRUE

  lengths <- ORFik:::optimizedTranscriptLengths(df)
  length <- lengths[tx_name == id]
  # Custom UTR lengths for Sac cer (yeast)
  if (organism(df) == "Saccharomyces cerevisiae")
    length$utr5_len <- 650

  start <- 1
  end <- tx_width

  if (length$cds_len > 0 & region_type %in% c("mrna", "leader+cds")) {
    start <- start + length$utr5_len
    end <- start + length$cds_len
  }
  grl <- GRangesList(GRanges("1", IRanges(start, end)))
  names(grl) <- id

  ranges <- unlistGrl(grl)
  ranges <- c(GRanges("1", IRanges(1, tx_width)), ranges)
  dt <- geneBoxFromRanges(ranges, tx_width,
                          cols = c("#FFFFFF", c("#F8766D","#00BA38","#619CFF")[start(ranges[-1]) %% 3 + 1]))[[1]]
  return(geneModelPanelPlot(dt))
}
