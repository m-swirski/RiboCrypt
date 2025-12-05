singlePlot_select_plot_type <- function(profile, withFrames, frame_colors, colors,
                                        lines, type, lib_index, line_size = 0.2) {
  count <- NULL # Avoid data.table warning
  profile_plot <- ggplot(profile)
  if (withFrames && type != "heatmap") frame_colors <- frame_color_themes(frame_colors)

  if (length(lines) > 0) profile_plot <- profile_plot +
    geom_vline(xintercept = lines, col = names(lines), linetype = 4,
               alpha = 0.2, size = line_size)
  profile_plot <- profile_plot + geom_hline(yintercept = 0)
  if (!withFrames) {
    profile_plot <- profile_plot +
      geom_area(aes(y = count, x = position), fill = colors, col = "black",
                size = 0.1, alpha = 0.8, position = "identity")
  } else if (type == "lines") {
    profile_plot <- profile_plot +
      geom_line(aes(y = count, x = position, color = frame), linewidth = 0.5) +
      scale_color_manual(values = frame_colors)
  } else if (type == "stacks") {
    profile_plot <- profile_plot +
      geom_area(aes(y = count, x = position, fill = frame), size = 0.1,
                alpha = 0.8, col = "black") +
      scale_fill_manual(values = frame_colors)
  } else if (type == "columns") {
    profile_plot <- profile_plot +
      geom_col(aes(y = count, x = position, fill = frame)) +
      scale_fill_manual(values = frame_colors)
  } else if (type == "area") {
    profile_plot <- profile_plot +
      geom_area(aes(y = count, x = position, fill = frame), size = 0.1,
                alpha = 0.8, col = "black", position = 'identity') +
      geom_point(aes(x = position, y = count, group = frame),
                 fill = NA,  # Make points invisible (no fill)
                 color = NA, size = 0.1, position = 'identity',
                 inherit.aes = FALSE) +
      scale_fill_manual(values = frame_colors)
  } else if (type == "heatmap") {
    hm_colors <- c("white", "yellow1","yellow2", "yellow3",
                   "lightblue", "blue", "navy")
    pro <- copy(profile)
    pro[, count := log2(count + 1)]
    profile_plot <- ggplot(pro, aes(y = lib_index, x = position, fill = count)) +
      geom_tile() +
      scale_fill_gradientn(colours = hm_colors,
                           limits = c(min(pro$count), max(pro$count)))
  }
  return(profile_plot)
}

singlePlot_add_theme <- function(profile_plot, ylabels, type,
                                 flip_ylabel = type == "heatmap", total_libs,
                                 ylabels_full_name = ylabels, as_plotly = TRUE) {
  y_text_size <- max(22 - total_libs * 3, 2)

  profile_plot <- profile_plot +
    ylab(ylabels) +
    theme(axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.y = element_text(size = y_text_size, face = "bold", angle = 0, vjust = 0.5),
          axis.text.y = element_text(size = y_text_size - 2),
          plot.margin = unit(c(0,0,0,0), "pt"),
          legend.position = "none") +
    scale_x_continuous(expand = c(0,0))

  if (type == "heatmap" | total_libs > 5) {
    profile_plot <- profile_plot +
      theme(axis.title.y = element_blank())
    if (type == "heatmap") {
      profile_plot <- profile_plot +
        theme(axis.ticks.y = element_blank(),
              axis.text.y = element_blank())
    }

    # browser()

    profile_plot <- automateTicksRNA(profile_plot, as_plotly, type != "heatmap",
                                     restrict_ticks = total_libs > 3)
    if (flip_ylabel | total_libs > 5) {
      y_text_size <- ifelse(total_libs < 30, 15, ifelse(total_libs < 50, 10,
                                                        ifelse(total_libs < 60, 7, 5)))
      # browser()
      profile_plot <- profile_plot %>%
        add_annotations(text = ylabels, x = "", y = 0.5,
                        yref = "paper", xref = "paper", showarrow = FALSE,
                        font = list(size = y_text_size),
                        hovertext = ylabels_full_name)
    }
    return(profile_plot)
  }
  return(automateTicksRNA(profile_plot, as_plotly, restrict_ticks = total_libs > 3))
}
