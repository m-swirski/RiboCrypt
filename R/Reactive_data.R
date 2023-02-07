heatmap_data <- function(mainPlotControls, tx, length_table) {
  message("-- Region: ", mainPlotControls()$region)
  if (length(mainPlotControls()$cds_display) > 0) {
    print("This is a mRNA")
    print(class(mainPlotControls()$reads[[1]]))
    # Pick start or stop region
    point <- observed_cds_point(mainPlotControls)
    windows <- startRegion(point, tx()[names(point)], TRUE,
                           upstream = mainPlotControls()$extendLeaders,
                           downstream = mainPlotControls()$extendTrailers - 1)
    length_table_sub <- length_table()[tx_name %in% names(point),]
    if (mainPlotControls()$region == "Start codon") {
      windows <- extend_needed(windows, length_table_sub$utr5_len,
                               mainPlotControls()$extendLeaders, "up")
      windows <- extend_needed(windows, length_table_sub$cds_len,
                               mainPlotControls()$extendTrailers - 1, "down")
    } else {
      windows <- extend_needed(windows, length_table_sub$cds_len,
                               mainPlotControls()$extendLeaders, "up")
      windows <- extend_needed(windows, length_table_sub$utr3_len,
                               mainPlotControls()$extendTrailers  - 1, "down")
    }
    window_lengths <- widthPerGroup(windows, FALSE)
    if (length(unique(window_lengths)) > 1) {
      warning("Some genes hit chromosome boundary, removing them.")
      windows <- windows[window_lengths == max(window_lengths)]
    }
    time_before <- Sys.time()
    # browser()

    dt <- windowPerReadLength(point, tx(),
                              reads = mainPlotControls()$reads[[1]],
                              pShifted = FALSE, upstream = mainPlotControls()$extendLeaders,
                              downstream = mainPlotControls()$extendTrailers - 1,
                              scoring = mainPlotControls()$normalization,
                              acceptedLengths = seq(mainPlotControls()$readlength_min, mainPlotControls()$readlength_max),
                              drop.zero.dt = TRUE, append.zeroes = TRUE,
                              windows = windows)
    if (!mainPlotControls()$p_shifted){

      sdt <- mainPlotControls()$shift_table
      if (nrow(sdt) > 0) {
        colnames(sdt)[1] <- "readlength"
        dt[, position := position + sdt[readlength == fraction]$offsets_start, by = fraction]
        dt <- dt[position %between% c(- mainPlotControls()$extendLeaders + max(abs(sdt$offsets_start)),
                                      mainPlotControls()$extendTrailers - 1 - max(abs(sdt$offsets_start)))]
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
