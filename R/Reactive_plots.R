reactive_heatmap_plot <- function(mainPlotControls, coverage) {
  {
    message("-- Plotting heatmap")
    pos <- ifelse(mainPlotControls()$region == "Start codon",
                  "Start Site", "Stop Site")
    main_plot <- coverageHeatMap(coverage(), scoring = mainPlotControls()$normalization,
                                 legendPos = "bottom",
                                 xlab = paste("Position relative to", pos)) +
      theme(axis.title = element_text(size = 18),
            axis.text = element_text(size = 12))
    plot_list <- if (mainPlotControls()$summary_track) {
      heights <- c(0.2,0.8)
      list(pSitePlot(coverage(), forHeatmap = TRUE), main_plot)
    } else {
      heights <- 1
      list(main_plot)
    }
    return(subplot(plot_list, nrows = length(plot_list), heights = heights, shareX = TRUE, titleY = TRUE) %>%
             plotly::config(toImageButtonOptions= list(format = "svg")))
  }
}
