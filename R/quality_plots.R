periodicity_plot <- function(dt, fft = FALSE) {
 score <- NULL # Avoid biocCheck error
 outplot <- ggplot(dt) +
    geom_col(aes(x = fraction, y = score, fill = frame), position = "dodge", width = .5) +
    scale_x_continuous(breaks=unique(dt$fraction)) +
    theme_bw() +
    theme(axis.text.y = element_blank(),
          axis.ticks = element_blank() ) +
   ylab("Frame score per read length")

  outplot <- ggplotly(outplot)
  if (fft) {
    rft <- dt_fft(dt)
    rft_p <-   lapply(unique(dt$fraction),
                      function(rl) ggplot(rft[read_length == rl & periods <= 5]) +
                        geom_line(aes(x = periods, y = amplitude)) +
                        theme_void()+
                        ylab("Amplitude"))
    rft_p <- subplot(rft_p, margin = 0)
    outplot <- subplot(outplot, rft_p, nrows = 2, margin = 0)
  }
  return(outplot)
}

