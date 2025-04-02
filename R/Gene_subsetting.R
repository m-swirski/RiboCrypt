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

extend_all_to <- function(point, tx, upstream, downstream) {
  windows <- startRegion(point, tx[names(point)], TRUE,
                         upstream,
                         downstream)

  TSS_to_anchor_dist <- pmapToTranscriptF(windows, tx[names(windows)])
  max_length_upstream <- unlist(end(TSS_to_anchor_dist), use.names = FALSE)
  max_length_downstream <- widthPerGroup(tx[names(windows)], keep.names = FALSE) - max_length_upstream

  windows <- extend_needed(windows, max_length_upstream, upstream, "up")
  windows <- extend_needed(windows, max_length_downstream, downstream, "down")

  window_lengths <- widthPerGroup(windows, FALSE)
  if (length(unique(window_lengths)) > 1) {
    warning("Some genes hit chromosome boundary, removing them.")
    table <- table(window_lengths)
    selected_length <- as.integer(names(which.max(table)))
    windows <- windows[window_lengths == selected_length]
  }
  return(windows)
}
