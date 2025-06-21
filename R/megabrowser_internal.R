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
                                  color_mult = 3, plotType = "plotly") {
  colors <- if (color_theme == "default (White-Blue)") {
    c("white", "lightblue", rep("blue", 4 + color_mult), "navy", "black")
  } else if (color_theme == "Matrix (black,green,red)") {
    c("#000000", "#2CFA1F", "yellow2", rep("#FF2400", color_mult))
  } else stop("Invalid color theme!")
  cat("Creating metabrowser heatmap\n")

  if (plotType == "plotly") {
    mat <- t(table)

    km <- kmeans(mat, centers = clusters)
    row_clusters <- seq_along(km$cluster)
    row_clusters <- split(row_clusters, km$cluster)

    if (clusters == 1) {
      row_clusters <- unlist(row_clusters, use.names = FALSE)
    }

    plot <- plotly::plot_ly(
      z = mat[unlist(row_clusters, use.names = FALSE),],
      colors = colors,
      type = "heatmapgl"
    )
    attr(plot, "row_order") <- row_clusters
  } else {
    plot <- ComplexHeatmap::Heatmap(t(table), show_row_dend = FALSE,
                            cluster_columns = FALSE,
                            cluster_rows = FALSE,
                            use_raster = TRUE,  raster_quality = 5,
                            km = clusters,
                            col =  colors, show_row_names = FALSE,
                            show_heatmap_legend = FALSE)
  }

  return(plot)

}

row_order <- function(x) {
  if (is(x, "Heatmap")) {
    ComplexHeatmap::row_order(x)
  } else attr(x, "row_order")
}

#' Full plot for allsamples browser
#'
#' @param m data.table of coverage per sample (wide format)
#' @param heatmap ComplexHeatmap object of plot from 'm'
#' @param id id of transcript
#' @param df ORFik experiment
#' @param summary logical, default TRUE (add top plot)
#' @param annotation logical, default TRUE (add bottom annotation track)
#' @param region_type character, "what is the coverage region?" Usually full mrna:
#' "mrna" or "leader+cds".
#' @param plotType = "plotly",
#' @param tx_annotation a GRangesList of tx annotation
#' @param display_region a GRangesList of display region
#' @param cds_annotation a GRangesList of cds annotation
#' @param viewMode character, "tx" or "genomic"
#' @param collapse_intron_flank integer, if TRUE and viewMode genomic, collapse
#' introns to this max size.
#' @param rel_heights numeric < 1, default: c(0.2, 0.75, 0.05).
#' Relative heights, sum to 1 and must be length 3.
#' @return a cowplot grub
#' @importFrom cowplot plot_grid
get_meta_browser_plot_full <- function(m, heatmap, id, df,
                                       summary = TRUE, annotation = TRUE,
                                       region_type,
                                       plotType = "plotly",
                                       tx_annotation, display_region,
                                       cds_annotation, viewMode,
                                       collapse_intron_flank,
                                       rel_heights = c(0.2, 0.75, 0.05)) {
  stopifnot(is.numeric(rel_heights) && length(rel_heights) == 3)
  gene_model_panel <- summary_plot <- NULL
  to_use_logicals <- c(summary, TRUE, annotation)

  if (!summary & !annotation) {
    gene_model_panel <- ggplot() + theme_classic()
    to_use_logicals[3] <- TRUE
  }

  if (summary) {
    summary_plot <- summary_track_allsamples(m)
  }

  if (annotation) {
    gene_model_panel <- annotation_track_allsamples(df, id, display_region, cds_annotation,
                                                    tx_annotation, viewMode, collapse_intron_flank)
  }
  if (plotType == "plotly") {
    gene_model_panel <- ggplotly(gene_model_panel, tooltip = "")
    final_plot <- heatmap
  } else {
    grob <- grid::grid.grabExpr(draw(heatmap))
    final_plot <- cowplot::plot_grid(plotlist = list(summary_plot, grob, gene_model_panel)[to_use_logicals],
                                     ncol = 1, rel_heights = rel_heights[to_use_logicals])
  }
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

annotation_track_allsamples <- function(df, id, display_range, annotation,
                                        tx_annotation, viewMode, collapse_intron_flank) {
  gene_model_panel <- createGeneModelPanel(display_range, annotation, frame = 1,
                                           tx_annotation = tx_annotation,
                                           custom_regions = NULL,
                                           viewMode = viewMode, collapse_intron_flank)
  return(geneModelPanelPlot(gene_model_panel[[1]]))
}
