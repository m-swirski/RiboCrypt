
timer_done_nice_print <- function(message = "", time_before) {
  cat(message); print(round(Sys.time() - time_before, 2))
}
