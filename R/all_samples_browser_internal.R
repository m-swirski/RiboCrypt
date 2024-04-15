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

get_meta_browser_plot_full <- function(m, heatmap, id, df,
                                       summary = TRUE, annotation = TRUE,
                                       rel_heights = c(0.2, 0.75, 0.05)) {
  stopifnot(is.numeric(rel_heights) && length(rel_heights) == 3)
  gene_model_panel <- summary_plot <- NULL
  if (!summary & !annotation) return(heatmap)
  # browser()
  if (summary) {
    summary_track_type <- "area"
    summary_profile <- data.table(count = rowSums(m))
    summary_profile[, `:=`(position = seq.int(.N)) ]
    summary_profile[, `:=`(frame = factor((position-1) %% 3)) ]
    summary_plot <- RiboCrypt:::createSinglePlot(summary_profile, TRUE, 1, "",
                                                 FALSE, lines = NULL,
                                                 type = summary_track_type,
                                                 flip_ylabel = FALSE, as_plotly = FALSE) +
      theme(axis.text.y=element_blank(),  axis.ticks.y=element_blank()) ; summary_plot
  }

  if (annotation) {
    lines <- NULL
    withFrames <- TRUE
    colors <- TRUE

    lengths <- ORFik:::optimizedTranscriptLengths(df)
    length <- lengths[tx_name == id]
    tx_width <- length$tx_len


    start <- 1
    end <- tx_width
    if (length$cds_len > 0) {
      start <- start + length$utr5_len
      end <- start + length$cds_len
    }
    grl <- GRangesList(GRanges("1", IRanges(start, end)))
    names(grl) <- id

    ranges <- unlistGrl(grl)
    ranges <- c(GRanges("1", IRanges(1, tx_width)), ranges)
    dt <- geneBoxFromRanges(ranges, tx_width,
                            cols = c("#FFFFFF", c("#F8766D","#00BA38","#619CFF")[start(ranges[-1]) %% 3 + 1]))[[1]]
    gene_model_panel <- geneModelPanelPlot(dt)
  }

  grob <- grid::grid.grabExpr(draw(heatmap))
  to_use_logicals <- c(summary, TRUE, annotation)
  final_plot <- cowplot::plot_grid(plotlist = list(summary_plot, grob, gene_model_panel)[to_use_logicals],
                                   ncol = 1, rel_heights = c(0.2, 0.75, 0.05)[to_use_logicals])
  return(final_plot)
}
