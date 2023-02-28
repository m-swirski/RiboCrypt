extend_needed <- function(windows, length, wanted_length, direction = "up") {
  dif <- length - wanted_length
  big_enough <- dif >= 0
  if (!all(big_enough)) {
    if (direction == "up") {
      windows[!big_enough] <- extendLeaders(windows[!big_enough],
                                            extension = -dif[!big_enough])
    } else {
      windows[!big_enough] <- extendTrailers(windows[!big_enough],
                                             extension = -dif[!big_enough])
    }

  }
  return(windows)
}

extend_all_to <- function(point, tx, length_table, mainPlotControls) {
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
    table <- table(window_lengths)
    selected_length <- as.integer(names(which.max(table)))
    windows <- windows[window_lengths == selected_length]
  }
  return(windows)
}
