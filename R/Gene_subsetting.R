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
