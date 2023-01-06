createSinglePlot <- function(profile, withFrames, colors, ylabels, lines, type = "lines"){
  profile_plot <- singlePlot_select_plot_type(profile, withFrames, colors,
                                              lines, type)
  return(singlePlot_add_theme(profile_plot, ylabels, type))
}




#profiles <- mapply(function(x,y,z) getProfileAnimate(display_range, x, y, z), reads, withFrames, kmers,  SIMPLIFY = FALSE)

getPlotAnimate <- function(profile, withFrames, colors, ylabels, lines){
  count <- NULL # Avoid data.table warning
  profile_plot <- ggplot(profile) +
    ylab(ylabels) +
    theme(axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          plot.margin = unit(c(0,0,0,0), "pt")) +
    scale_x_continuous(expand = c(0,0))

  if (withFrames) {
    if (length(lines) > 0) profile_plot <- profile_plot + geom_vline(xintercept = lines, col = names(lines), linetype = 4)
    profile_plot <- profile_plot +
      geom_line(aes(y = count, x = position, color = frame, frame = file), size = 0.75) +
      theme(legend.position = "none") +
      xlab("position [nt]")
    profile_plot <- automateTicks(profile_plot)

  } else {
    if (length(lines) > 0) profile_plot <- profile_plot + geom_vline(xintercept = lines, col = names(lines), linetype = 4)
    profile_plot <- profile_plot +
      geom_area(aes(y = count, x = position, frame = file), fill = colors, position = "identity")
    profile_plot <- automateTicksRNA(profile_plot)
  }
  return(profile_plot)
}
