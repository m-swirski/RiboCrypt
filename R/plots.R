createSinglePlot <- function(profile, withFrames, frame_colors, colors, ylabels,
                             ylabels_full_name = ylabels, lines, type = "lines",
                             lib_index, total_libs = 1,
                             flip_ylabel = type == "heatmap",
                             as_plotly = TRUE,
                             templates = NULL) {
  profile_plot <- singlePlot_select_plot_type(
    profile, withFrames, frame_colors, colors,
    lines, type, lib_index, templates = templates
  )
  singlePlot_add_theme(
    profile_plot, ylabels, type, flip_ylabel,
    total_libs, ylabels_full_name, as_plotly
  )
}

getPlotAnimate <- function(profile, withFrames, colors, frame_colors,
                           ylabels, lines, lines_size = 0.1) {
  profile <- data.table::as.data.table(profile)

  if (withFrames) {
    frame_levels <- track_frame_levels(profile)
    frame_colors <- track_frame_colors(frame_levels, frame_colors)
    plot <- suppressWarnings(plotly::plot_ly(
      profile,
      x = ~position,
      y = ~count,
      frame = ~file,
      split = ~frame,
      type = "scatter",
      mode = "lines",
      colors = unname(frame_colors),
      color = ~frame,
      line = list(simplify = FALSE),
      hovertemplate = "position: %{x:.0f}<br>count: %{y}<br>frame: %{fullData.name}<extra></extra>"
    ))
  } else {
    plot <- suppressWarnings(plotly::plot_ly(
      profile,
      x = ~position,
      y = ~count,
      frame = ~file,
      type = "scatter",
      mode = "lines",
      line = list(color = "black", width = 0.5, simplify = FALSE),
      fill = "tozeroy",
      fillcolor = grDevices::adjustcolor(
        resolve_track_fill_color(colors),
        alpha.f = track_area_fill_alpha()
      ),
      hovertemplate = "position: %{x:.0f}<br>count: %{y}<extra></extra>",
      showlegend = FALSE
    ))
  }

  y_max <- max(profile$count, na.rm = TRUE)
  if (!is.finite(y_max) || y_max <= 0) y_max <- 1
  plot <- track_guides(
    plot,
    x_range = range(profile$position, na.rm = TRUE),
    y_max = y_max,
    lines = lines,
    add_zero = FALSE,
    line_size = lines_size,
    line_alpha = 0.5
  )

  plot <- plot %>%
    plotly::layout(
      xaxis = list(
        autorange = FALSE,
        title = list(text = ""),
        showticklabels = FALSE,
        ticks = ""
      ),
      yaxis = list(
        autorange = TRUE,
        title = list(text = ylabels),
        zeroline = FALSE
      ),
      margin = list(t = 0, r = 0, b = 0, l = 0, pad = 0)
    )
  suppressWarnings(plotly::animation_opts(plot, frame = 80, transition = 0, redraw = FALSE))
}

make_summary_track <- function(profiles, plots, withFrames, frame_colors, colors,
                               lines, summary_track_type, nplots, templates = NULL) {
  count <- NULL # avoid bioccheck error
  summary_profile <- data.table::rbindlist(profiles)
  summary_profile <- summary_profile[, .(count = sum(count)), by = position]
  if (!is.null(profiles[[1]]$frame)) summary_profile[, frame := profiles[[1]]$frame]
  summary_plot <- createSinglePlot(
    summary_profile, all(withFrames), frame_colors, colors[1], "summary",
    FALSE, lines, type = summary_track_type, flip_ylabel = FALSE,
    lib_index = nplots, total_libs = nplots, templates = templates
  )
  plots[[nplots]] <- summary_plot
  rev(plots)
}
