createSinglePlot <- function(profile, withFrames, frame_colors, colors, ylabels,
                             ylabels_full_name = ylabels, lines, type = "lines",
                             lib_index, total_libs = 1,
                             flip_ylabel = type == "heatmap",
                             as_plotly = TRUE){
  profile_plot <- singlePlot_select_plot_type(profile, withFrames, frame_colors, colors,
                                              lines, type, lib_index)
  return(singlePlot_add_theme(profile_plot, ylabels, type, flip_ylabel,
                              total_libs, ylabels_full_name, as_plotly))
}




#profiles <- mapply(function(x,y,z) getProfileAnimate(display_range, x, y, z), reads, withFrames, kmers,  SIMPLIFY = FALSE)

getPlotAnimate <- function(profile, withFrames, colors, frame_colors,
                           ylabels, lines, lines_size = 0.1){
  count <- NULL # Avoid data.table warning
  profile_plot <- ggplot(profile) +
    ylab(ylabels) +
    theme(axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          plot.margin = unit(c(0,0,0,0), "pt")) +
    scale_x_continuous(expand = c(0,0))
  if (withFrames) colors <- frame_color_themes(frame_colors)

  if (length(lines) > 0) profile_plot <- profile_plot +
    geom_vline(xintercept = lines, col = names(lines), linetype = 4, alpha = 0.5, size = lines_size)
  if (withFrames) {
    profile_plot <- profile_plot +
      suppressWarnings(geom_line(aes(y = count, x = position, color = frame, frame = file), size = 0.75)) +
      scale_color_manual(values = colors) +
      theme(legend.position = "none") +
      xlab("position [nt]")
    profile_plot <- automateTicks(profile_plot)

  } else {
    profile_plot <- profile_plot +
      suppressWarnings(geom_area(aes(y = count, x = position, frame = file), fill = colors, position = "identity"))
    profile_plot <- automateTicksRNA(profile_plot)
  }
  return(profile_plot)
}

make_summary_track <- function(profiles, plots, withFrames, frame_colors, colors,
                               lines, summary_track_type, nplots) {
  count <- NULL # avoid bioccheck error
  summary_profile <- rbindlist(profiles)
  summary_profile <- summary_profile[,.(count = sum(count)), by = position]
  if (!is.null(profiles[[1]]$frame)) summary_profile[, frame := profiles[[1]]$frame]
  summary_plot <- createSinglePlot(summary_profile, all(withFrames), frame_colors, colors[1], "summary",
                                   FALSE, lines,
                                   type = summary_track_type,
                                   flip_ylabel = FALSE)
  plots[[nplots]] <- summary_plot
  plots <- rev(plots)
  return(plots)
}
