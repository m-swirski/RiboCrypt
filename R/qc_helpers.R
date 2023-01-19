
dt_fft <- function(dt) {
  read_lengths <- unique(dt$fraction)
  fft_dt <- data.table()
  for (i in read_lengths) {
    spec <- spec.pgram(x = dt[fraction == i, ]$score, plot = FALSE)
    fft_dt <- rbindlist(list(fft_dt, data.table(read_length = i, 
                                                amplitude = spec$spec, periods = 1/spec$freq)))
  }
  return(fft_dt)
}
