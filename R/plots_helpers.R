singlePlot_select_plot_type <- function(profile, withFrames, colors,
                                        lines, type, line_size = 0.2) {
  count <- NULL # Avoid data.table warning
  profile_plot <- ggplot(profile)
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
      geom_line(aes(y = count, x = position, color = frame), size = 0.5)
  } else if (type == "stacks") {
    profile_plot <- profile_plot +
      geom_area(aes(y = count, x = position, fill = frame), size = 0.1,
                alpha = 0.8, col = "black")
  } else if (type == "columns") {
    profile_plot <- profile_plot +
      geom_col(aes(y = count, x = position, fill = frame))
  } else if (type == "area") {
    profile_plot <- profile_plot +
      geom_area(aes(y = count, x = position, fill = frame), size = 0.1,
                alpha = 0.8, col = "black", position = 'identity') +
      geom_point(aes(x = position, y = count, group = frame),
                 fill = NA,  # Make points invisible (no fill)
                 color = NA, size = 0.1, position = 'identity',
                 inherit.aes = FALSE)
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

singlePlot_add_theme <- function(profile_plot, ylabels, type,
                                 flip_ylabel = type == "heatmap", total_libs,
                                 ylabels_full_name = ylabels, as_plotly = TRUE,
                                 y_autorange = FALSE, y_nticks = 3) {
  y_text_size <- ifelse(total_libs == 1, 8, ifelse(total_libs < 30, 7, 6))
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
    profile_plot <- if (type == "heatmap") {
      automateTicksRNA(profile_plot, as_plotly)
    } else automateTicksRNA(profile_plot, as_plotly, y_autorange, y_nticks)
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
  return(automateTicksRNA(profile_plot, as_plotly, y_autorange, y_nticks))
}
