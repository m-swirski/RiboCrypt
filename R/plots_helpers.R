resolve_track_color <- function(color) {
  if (length(color) == 0 || is.null(color) || is.na(color[[1]])) return("#4C78A8")
  color <- color[[1]]
  if (is.numeric(color)) {
    palette_colors <- grDevices::palette()
    index <- ((as.integer(color) - 1L) %% length(palette_colors)) + 1L
    return(palette_colors[[index]])
  }
  as.character(color)
}

track_frame_levels <- function(profile) {
  frame <- profile$frame
  if (is.factor(frame)) {
    levels(frame)[levels(frame) %in% as.character(unique(frame))]
  } else {
    sort(unique(as.character(frame)))
  }
}

track_frame_colors <- function(frame_levels, frame_colors, with_alpha = FALSE) {
  colors <- frame_color_themes(frame_colors, with_alpha)
  stats::setNames(colors[seq_along(frame_levels)], frame_levels)
}

resolve_track_fill_color <- function(color) {
  alpha_normalize_colors(resolve_track_color(color))[[1]]
}

track_hover_template <- function(include_frame = FALSE, type = "scatter") {
  base <- "position: %{x:.0f}<br>count: %{y}"
  if (type == "heatmap") {
    base <- "position: %{x:.0f}<br>count: %{z}"
  }
  if (include_frame) base <- paste0(base, "<br>frame: %{fullData.name}")
  paste0(base, "<extra></extra>")
}

track_area_fill_alpha <- function() 0.6

track_guides <- function(plot, x_range, y_max, lines = numeric(), add_zero = TRUE,
                         line_size = 0.2, line_alpha = 0.2) {
  if (length(lines) > 0) {
    for (i in seq_along(lines)) {
      plot <- plotly::add_trace(
        plot,
        x = c(unname(lines[[i]]), unname(lines[[i]])),
        y = c(0, y_max),
        type = "scatter",
        mode = "lines",
        line = list(
          color = names(lines)[[i]],
          width = line_size,
          dash = "dot"
        ),
        opacity = line_alpha,
        hoverinfo = "skip",
        showlegend = FALSE,
        inherit = FALSE
      )
    }
  }
  if (isTRUE(add_zero)) {
    plot <- plotly::add_trace(
      plot,
      x = x_range,
      y = c(0, 0),
      type = "scatter",
      mode = "lines",
      line = list(color = "rgba(0,0,0,0.5)", width = 0.5),
      hoverinfo = "skip",
      showlegend = FALSE,
      inherit = FALSE
    )
  }
  plot
}

singlePlot_select_plot_type <- function(profile, withFrames, frame_colors, colors,
                                        lines, type, lib_index, line_size = 0.2) {
  count <- NULL # Avoid data.table warning
  profile <- data.table::as.data.table(profile)
  frame_theme <- frame_colors

  if (withFrames && type != "heatmap") {
    frame_levels <- track_frame_levels(profile)
    frame_colors <- track_frame_colors(frame_levels, frame_theme)
  } else {
    frame_levels <- character()
  }

  if (type == "heatmap") {
    hm_colors <- c("white", "yellow1", "yellow2", "yellow3",
                   "lightblue", "blue", "navy")
    pro <- data.table::copy(profile)
    pro[, count := log2(count + 1)]

    return(plotly::plot_ly(
      x = pro$position,
      y = lib_index,
      z = matrix(pro$count, nrow = 1),
      type = "heatmap",
      colors = hm_colors,
      showscale = FALSE,
      hovertemplate = track_hover_template(type = "heatmap"),
      showlegend = FALSE
    ))
  }

  plot <- plotly::plot_ly()
  guide_y_max <- max(profile$count, na.rm = TRUE)
  x_range <- range(profile$position, na.rm = TRUE)
  if (!is.finite(guide_y_max) || guide_y_max <= 0) guide_y_max <- 1

  if (!withFrames) {
    fill_color <- resolve_track_fill_color(colors)
    plot <- plotly::add_trace(
      plot,
      data = profile,
      x = ~position,
      y = ~count,
      type = "scatter",
      mode = "lines",
      line = list(color = "black", width = 0.5),
      fill = "tozeroy",
      fillcolor = grDevices::adjustcolor(fill_color, alpha.f = track_area_fill_alpha()),
      hovertemplate = track_hover_template(),
      showlegend = FALSE
    )
  } else if (type == "lines") {
    for (frame_name in frame_levels) {
      frame_dt <- profile[as.character(frame) == frame_name]
      plot <- plotly::add_trace(
        plot,
        data = frame_dt,
        x = ~position,
        y = ~count,
        type = "scatter",
        mode = "lines",
        name = frame_name,
        legendgroup = frame_name,
        line = list(color = frame_colors[[frame_name]], width = 0.8),
        hovertemplate = track_hover_template(include_frame = TRUE),
        showlegend = TRUE
      )
    }
  } else if (type == "columns") {
    for (frame_name in frame_levels) {
      frame_dt <- profile[as.character(frame) == frame_name]
      plot <- plotly::add_trace(
        plot,
        data = frame_dt,
        x = ~position,
        y = ~count,
        type = "bar",
        name = frame_name,
        legendgroup = frame_name,
        marker = list(color = frame_colors[[frame_name]]),
        hovertemplate = track_hover_template(include_frame = TRUE),
        showlegend = TRUE
      )
    }
  } else if (type == "area") {
    frame_fill_colors <- track_frame_colors(frame_levels, frame_theme, with_alpha = TRUE)
    for (frame_name in frame_levels) {
      frame_dt <- profile[as.character(frame) == frame_name]
      plot <- plotly::add_trace(
        plot,
        data = frame_dt,
        x = ~position,
        y = ~count,
        type = "scatter",
        mode = "lines",
        name = frame_name,
        legendgroup = frame_name,
        line = list(color = "black", width = 0.4),
        fill = "tozeroy",
        fillcolor = grDevices::adjustcolor(
          frame_fill_colors[[frame_name]],
          alpha.f = track_area_fill_alpha()
        ),
        hovertemplate = track_hover_template(include_frame = TRUE),
        showlegend = TRUE
      )
    }
  } else if (type == "stacks") {
    for (frame_name in frame_levels) {
      frame_dt <- profile[as.character(frame) == frame_name]
      plot <- plotly::add_trace(
        plot,
        data = frame_dt,
        x = ~position,
        y = ~count,
        type = "scatter",
        mode = "lines",
        name = frame_name,
        legendgroup = frame_name,
        stackgroup = "coverage",
        line = list(color = "black", width = 0.4),
        fillcolor = grDevices::adjustcolor(frame_colors[[frame_name]], alpha.f = 0.8),
        hovertemplate = track_hover_template(include_frame = TRUE),
        showlegend = TRUE
      )
    }
    guide_y_max <- max(profile[, .(count = sum(count)), by = position]$count, na.rm = TRUE)
  }

  track_guides(plot, x_range, guide_y_max, lines, add_zero = TRUE, line_size = line_size)
}

singlePlot_add_theme <- function(profile_plot, ylabels, type,
                                 flip_ylabel = type == "heatmap", total_libs,
                                 ylabels_full_name = ylabels, as_plotly = TRUE) {
  y_text_size <- max(22 - total_libs * 3, 2)
  annotation_list <- list()

  yaxis <- list(
    autorange = TRUE,
    rangemode = if (type == "heatmap") "normal" else "tozero",
    zeroline = type != "heatmap",
    title = list(text = ylabels, font = list(size = y_text_size)),
    tickfont = list(size = 16)
  )

  xaxis <- list(
    autorange = FALSE,
    title = list(text = ""),
    showticklabels = FALSE,
    ticks = "",
    showgrid = FALSE,
    zeroline = FALSE
  )

  layout_args <- list(
    xaxis = xaxis,
    yaxis = yaxis,
    margin = list(t = 0, r = 0, b = 0, l = 0, pad = 0),
    legend = list(orientation = "h"),
    hovermode = "x unified"
  )

  if (type == "columns") {
    layout_args$barmode <- "stack"
    layout_args$bargap <- 0.12
  }

  if (type == "heatmap" || total_libs > 5) {
    layout_args$yaxis$title <- NULL
    if (type == "heatmap") {
      layout_args$yaxis$showticklabels <- FALSE
      layout_args$yaxis$ticks <- ""
    } else if (is.numeric(ylabels) || grepl("^[0-9]+$", ylabels)) {
      layout_args$yaxis$tickmode <- "array"
      layout_args$yaxis$tickvals <- 0
      layout_args$yaxis$ticktext <- ""
    }

    if (flip_ylabel || total_libs > 5) {
      y_text_size <- ifelse(total_libs < 30, 15, ifelse(total_libs < 50, 10,
                                                        ifelse(total_libs < 60, 7, 5)))
      annotation_list <- list(list(
        text = ylabels,
        x = 0,
        y = 0.5,
        xref = "paper",
        yref = "paper",
        xanchor = "right",
        yanchor = "middle",
        showarrow = FALSE,
        font = list(size = y_text_size),
        hovertext = ylabels_full_name
      ))
    }
  } else {
    layout_args$yaxis$nticks <- 3
    if (total_libs > 3) layout_args$yaxis$nticks <- 2
  }

  if (length(annotation_list) > 0) layout_args$annotations <- annotation_list

  do.call(plotly::layout, c(list(profile_plot), layout_args))
}
