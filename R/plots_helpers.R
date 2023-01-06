singlePlot_select_plot_type <- function(profile, withFrames, colors,
                                        lines, type, line_size = 0.2) {
  count <- NULL # Avoid data.table warning
  profile_plot <- ggplot(profile)
  if (length(lines) > 0) profile_plot <- profile_plot +
    geom_vline(xintercept = lines, col = names(lines), linetype = 4, size = line_size)
  profile_plot <- profile_plot + geom_hline(yintercept = 0)
  if (!withFrames) {
    profile_plot <- profile_plot +
      geom_area(aes(y = count, x = position), fill = colors, col = "black",
                size = 0.1, alpha = 0.8, position = "identity")
  } else if (type == "lines") {
    profile_plot <- profile_plot +
      geom_line(aes(y = count, x = position,color = frame), size = 0.5)
  } else if (type == "stacks") {
    profile_plot <- profile_plot +
      geom_area(aes(y = count, x = position,fill = frame), size = 0.1,
                alpha = 0.8, col = "black")
  } else if (type == "columns") {
    profile_plot <- profile_plot +
      geom_col(aes(y = count, x = position, fill = frame))
  } else if (type == "area") {
    profile_plot <- profile_plot +
      geom_area(aes(y = count, x = position,fill = frame), size = 0.1,
                alpha = 0.8, col = "black", position = 'identity')
  } else if (type == "heatmap") {
    hm_colors <- c("white", "yellow1","yellow2", "yellow3",
                   "lightblue", "blue", "navy")
    pro <- copy(profile)
    pro[, count := log2(count + 1)]
    profile_plot <- ggplot(pro, aes(y = colors, x = position, fill = count)) +
      geom_tile() +
      scale_fill_gradientn(colours = hm_colors,
                           limits = c(min(pro$count), max(pro$count)))
  }
  return(profile_plot)
}

singlePlot_add_theme <- function(profile_plot, ylabels, type) {
  if (type == "heatmap") {
    profile_plot <- profile_plot +
      theme(legend.position = "none") +
      ylab(ylabels) +
      theme(axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            axis.title.y = element_text(size = 8, face = "bold"),
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank()) +
      theme(plot.margin = unit(c(0,0,0,0), "pt"))+
      scale_x_continuous(expand = c(0,0))
  } else {
    profile_plot <- profile_plot +
      theme(legend.position = "none") +
      ylab(ylabels) +
      theme(axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            axis.title.y = element_text(size = 8, face = "bold"),
            axis.text.y = element_text(size = 6)) +
      theme(plot.margin = unit(c(0,0,0,0), "pt"))+
      scale_x_continuous(expand = c(0,0))
  }

  return(automateTicksRNA(profile_plot))
}