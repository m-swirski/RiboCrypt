createSinglePlot <- function(profile, withFrames, colors, ylabels,
                             ylabels_full_name = ylabels, lines, type = "lines",
                             total_libs = 1, flip_ylabel = type == "heatmap"){
  profile_plot <- singlePlot_select_plot_type(profile, withFrames, colors,
                                              lines, type)
  return(singlePlot_add_theme(profile_plot, ylabels, type, flip_ylabel,
                              total_libs, ylabels_full_name))
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

make_summary_track <- function(profiles, plots, withFrames, colors, lines, summary_track_type, nplots) {
  count <- NULL # avoid bioccheck error
  summary_profile <- rbindlist(profiles)
  summary_profile <- summary_profile[,.(count = sum(count)), by = position]
  summary_profile[, frame := profiles[[1]]$frame]
  summary_plot <- createSinglePlot(summary_profile, all(withFrames), colors[1], "summary",
                                   FALSE, lines,
                                   type = summary_track_type,
                                   flip_ylabel = FALSE)
  plots[[nplots]] <- summary_plot
  plots <- rev(plots)
  return(plots)
}
