heatmap_data <- function(mainPlotControls, tx, anchor_points) {
  message("-- Region (motif anchor): ", mainPlotControls()$region)
  if (length(mainPlotControls()$cds_display) > 0) {

    windows <- extend_all_to(anchor_points, tx(),
                             mainPlotControls()$extendLeaders,
                             mainPlotControls()$extendTrailers - 1)
    time_before <- Sys.time()
    dt <- windowPerReadLength(anchor_points, tx(),
                              reads = mainPlotControls()$reads[[1]],
                              pShifted = FALSE, upstream = mainPlotControls()$extendLeaders,
                              downstream = mainPlotControls()$extendTrailers - 1,
                              scoring = mainPlotControls()$normalization,
                              acceptedLengths = seq(mainPlotControls()$readlength_min, mainPlotControls()$readlength_max),
                              drop.zero.dt = TRUE, append.zeroes = TRUE,
                              windows = windows)
    if (!is.null(mainPlotControls()$p_shifted)) {
      if (!mainPlotControls()$p_shifted){
        sdt <- mainPlotControls()$shift_table
        if (nrow(sdt) > 0) {
          colnames(sdt)[1] <- "readlength"
          dt[, position := position + sdt[readlength == fraction]$offsets_start, by = fraction]
          dt <- dt[position %between% c(- mainPlotControls()$extendLeaders + max(abs(sdt$offsets_start)),
                                        mainPlotControls()$extendTrailers - 1 - max(abs(sdt$offsets_start)))]
        }
      }
    }
    print(paste("Number of rows in dt:", nrow(dt)))
    cat("Coverage calc: "); print(round(Sys.time() - time_before, 2))
    return(dt)
  } else {
    print("This is not a mRNA / valid mRNA")
    return(NULL)
  }
}
