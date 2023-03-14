heatmap_data <- function(mainPlotControls, tx, length_table) {
  message("-- Region: ", mainPlotControls()$region)
  if (length(mainPlotControls()$cds_display) > 0) {
    print(paste("Number of input ranges: ",
                length(mainPlotControls()$cds_display)))
    print(class(mainPlotControls()$reads[[1]]))
    # Pick start or stop region
    point <- observed_cds_point(mainPlotControls)
    windows <- extend_all_to(point, tx, length_table, mainPlotControls)
    time_before <- Sys.time()
    dt <- windowPerReadLength(point, tx(),
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
