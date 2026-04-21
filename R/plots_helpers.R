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

track_vertical_segments <- function(x, y, text = NULL) {
  xseg <- c(rbind(x, x, NA_real_))
  yseg <- c(rbind(rep(0, length(y)), y, NA_real_))
  textseg <- NULL
  if (!is.null(text)) {
    textseg <- c(rbind(text, text, NA_character_))
  }
  list(x = xseg, y = yseg, text = textseg)
}

ensure_columns_gl_hover <- function(plot) {
  for (i in seq_along(plot$x$data)) {
    trace <- plot$x$data[[i]]
    is_columns_gl_trace <- identical(trace$type, "scattergl") &&
      identical(trace$mode, "lines")
    if (!is_columns_gl_trace) next
    plot$x$data[[i]]$hoverinfo <- "text"
    plot$x$data[[i]]$hovertemplate <- "%{text}<extra></extra>"
  }
  plot
}

columns_zoom_switch_threshold <- function() 5000

covPanelWithFramesPlotlyTemplate <- function() {
  frame_names <- c("0", "1", "2")
  p <- plotly::plot_ly()
  for (frame_name in frame_names) {
    p <- plotly::add_trace(
      p,
      x = numeric(),
      y = numeric(),
      type = "scatter",
      mode = "lines",
      name = frame_name,
      legendgroup = frame_name,
      line = list(color = "#000000", width = 0.8, simplify = FALSE),
      hovertemplate = track_hover_template(include_frame = TRUE),
      showlegend = TRUE,
      inherit = FALSE
    )
  }
  plotly::plotly_build(p)
}

covPanelWithoutFramesPlotlyTemplate <- function() {
  plotly::plotly_build(
    plotly::plot_ly(
      x = numeric(),
      y = numeric(),
      type = "scatter",
      mode = "lines",
      line = list(color = "black", width = 0.5, simplify = FALSE),
      fill = "tozeroy",
      fillcolor = resolve_track_fill_color(NULL),
      hovertemplate = track_hover_template(),
      showlegend = FALSE
    )
  )
}

covPanelColumnsPlotlyTemplate <- function() {
  p <- plotly::plot_ly()
  for (frame_name in c("0", "1", "2")) {
    p <- plotly::add_trace(
      p,
      x = c(NA_real_),
      y = c(NA_real_),
      type = "scattergl",
      mode = "lines",
      name = frame_name,
      legendgroup = frame_name,
      line = list(color = "#000000", width = 6, simplify = FALSE),
      text = NA_character_,
      meta = list(rc_columns_switch = "gl_subset"),
      showlegend = TRUE,
      visible = "legendonly",
      inherit = FALSE
    )
  }
  for (frame_name in c("0", "1", "2")) {
    p <- plotly::add_trace(
      p,
      x = numeric(),
      y = numeric(),
      type = "scatter",
      mode = "lines",
      name = frame_name,
      legendgroup = frame_name,
      line = list(color = "#000000", width = 0.8, simplify = FALSE),
      hovertemplate = track_hover_template(include_frame = TRUE),
      meta = list(rc_columns_switch = "line"),
      showlegend = FALSE,
      visible = TRUE,
      inherit = FALSE
    )
  }
  ensure_columns_gl_hover(plotly::plotly_build(p))
}

covPanelColumnsGLPlotlyTemplate <- function(bar_px = 6) {
  p <- plotly::plot_ly()
  for (frame_name in c("0", "1", "2")) {
    p <- plotly::add_trace(
      p,
      x = c(NA_real_),
      y = c(NA_real_),
      type = "scattergl",
      mode = "lines",
      name = frame_name,
      legendgroup = frame_name,
      line = list(color = "#000000", width = bar_px, simplify = FALSE),
      text = NA_character_,
      showlegend = TRUE,
      inherit = FALSE
    )
  }
  ensure_columns_gl_hover(plotly::plotly_build(p))
}

covPanelAreaPlotlyTemplate <- function() {
  p <- plotly::plot_ly()
  for (frame_name in c("0", "1", "2")) {
    p <- plotly::add_trace(
      p,
      x = numeric(),
      y = numeric(),
      type = "scatter",
      mode = "lines",
      name = frame_name,
      legendgroup = frame_name,
      line = list(color = "black", width = 0.4, simplify = FALSE),
      fill = "tozeroy",
      fillcolor = "#00000099",
      hovertemplate = track_hover_template(include_frame = TRUE),
      showlegend = TRUE,
      inherit = FALSE
    )
  }
  plotly::plotly_build(p)
}

covPanelHeatmapPlotlyTemplate <- function() {
  hm_colors <- c("white", "yellow1", "yellow2", "yellow3",
                 "lightblue", "blue", "navy")
  plotly::plotly_build(
    plotly::plot_ly(
      x = 1,
      y = 1,
      z = matrix(0, nrow = 1),
      type = "heatmap",
      colors = hm_colors,
      showscale = FALSE,
      hovertemplate = track_hover_template(type = "heatmap"),
      showlegend = FALSE
    )
  )
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
        line = list(color = "rgba(0,0,0,0.5)", width = 0.5, simplify = FALSE),
        hoverinfo = "skip",
        showlegend = FALSE,
        inherit = FALSE
    )
  }
  plot
}

singlePlot_select_plot_type <- function(profile, withFrames, frame_colors, colors,
                                        lines, type, lib_index, line_size = 0.2,
                                        templates = NULL) {
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
    pro <- data.table::copy(profile)
    pro[, count := log2(count + 1)]
    if (inherits(templates$cov_panel_heatmap_plotly, "plotly")) {
      plot <- templates$cov_panel_heatmap_plotly
      plot$x$data[[1]]$x <- pro$position
      plot$x$data[[1]]$y <- lib_index
      plot$x$data[[1]]$z <- matrix(pro$count, nrow = 1)
      plot$x$data[[1]]$colorscale <- NULL
      plot$x$data[[1]]$zmin <- NULL
      plot$x$data[[1]]$zmax <- NULL
      plot$x$data[[1]]$zauto <- NULL
      plot$x$data[[1]]$autocolorscale <- NULL
      if (length(plot$x$attrs) >= 1) {
        plot$x$attrs[[1]]$x <- pro$position
        plot$x$attrs[[1]]$y <- lib_index
        plot$x$attrs[[1]]$z <- matrix(pro$count, nrow = 1)
      }
      return(plot)
    }

    hm_colors <- c("white", "yellow1", "yellow2", "yellow3",
                   "lightblue", "blue", "navy")
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
  columns_trace_group <- paste0("track_", lib_index)
  guide_y_max <- max(profile$count, na.rm = TRUE)
  x_range <- range(profile$position, na.rm = TRUE)
  if (!is.finite(guide_y_max) || guide_y_max <= 0) guide_y_max <- 1

  if (!withFrames) {
    if (type == "lines" && inherits(templates$cov_panel_without_frames_plotly, "plotly")) {
      plot <- templates$cov_panel_without_frames_plotly
      fill_color <- resolve_track_fill_color(colors)
      plot$x$data[[1]]$x <- profile$position
      plot$x$data[[1]]$y <- profile$count
      plot$x$data[[1]]$fillcolor <- grDevices::adjustcolor(
        fill_color,
        alpha.f = track_area_fill_alpha()
      )
    } else {
      fill_color <- resolve_track_fill_color(colors)
      plot <- plotly::add_trace(
        plot,
        data = profile,
        x = ~position,
        y = ~count,
        type = "scatter",
        mode = "lines",
        line = list(color = "black", width = 0.5, simplify = FALSE),
        fill = "tozeroy",
        fillcolor = grDevices::adjustcolor(fill_color, alpha.f = track_area_fill_alpha()),
        hovertemplate = track_hover_template(),
        showlegend = FALSE
      )
    }
  } else if (type == "lines") {
    if (inherits(templates$cov_panel_with_frames_plotly, "plotly")) {
      plot <- templates$cov_panel_with_frames_plotly
      for (i in seq_along(frame_levels)) {
        frame_name <- frame_levels[[i]]
        frame_dt <- profile[as.character(frame) == frame_name]
        plot$x$data[[i]]$x <- frame_dt$position
        plot$x$data[[i]]$y <- frame_dt$count
        plot$x$data[[i]]$name <- frame_name
        plot$x$data[[i]]$legendgroup <- frame_name
        plot$x$data[[i]]$line$color <- frame_colors[[frame_name]]
      }
      if (length(plot$x$data) > length(frame_levels)) {
        plot$x$data <- plot$x$data[seq_along(frame_levels)]
      }
    } else {
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
          line = list(color = frame_colors[[frame_name]], width = 0.8, simplify = FALSE),
          hovertemplate = track_hover_template(include_frame = TRUE),
          showlegend = TRUE
        )
      }
    }
  } else if (type == "columns") {
    if (inherits(templates$cov_panel_columns_plotly, "plotly")) {
      plot <- templates$cov_panel_columns_plotly
      is_gl_template <- identical(plot$x$data[[1]]$type, "scattergl")
      has_switch_lines <- length(plot$x$data) > length(frame_levels)
      for (i in seq_along(frame_levels)) {
        frame_name <- frame_levels[[i]]
        frame_dt <- profile[as.character(frame) == frame_name]
        plot$x$data[[i]]$name <- frame_name
        plot$x$data[[i]]$legendgroup <- frame_name
        if (is_gl_template) {
          plot$x$data[[i]]$line$color <- frame_colors[[frame_name]]
          plot$x$data[[i]]$meta <- c(
            plot$x$data[[i]]$meta,
            list(
              rc_columns_group = columns_trace_group,
              rc_columns_frame = frame_name
            )
          )
          if (has_switch_lines) {
            plot$x$data[[i]]$x <- c(NA_real_)
            plot$x$data[[i]]$y <- c(NA_real_)
            plot$x$data[[i]]$text <- NA_character_
            plot$x$data[[i]]$visible <- "legendonly"
            line_index <- i + length(frame_levels)
            plot$x$data[[line_index]]$x <- frame_dt$position
            plot$x$data[[line_index]]$y <- frame_dt$count
            plot$x$data[[line_index]]$line$color <- frame_colors[[frame_name]]
            plot$x$data[[line_index]]$meta <- c(
              plot$x$data[[line_index]]$meta,
              list(
                rc_columns_group = columns_trace_group,
                rc_columns_frame = frame_name
              )
            )
            plot$x$data[[line_index]]$visible <- TRUE
          } else {
            hover_text <- paste0(
              "position: ", frame_dt$position,
              "<br>count: ", frame_dt$count,
              "<br>frame: ", frame_name
            )
            seg <- track_vertical_segments(frame_dt$position, frame_dt$count, hover_text)
            plot$x$data[[i]]$x <- seg$x
            plot$x$data[[i]]$y <- seg$y
            plot$x$data[[i]]$text <- seg$text
          }
        } else {
          plot$x$data[[i]]$x <- frame_dt$position
          plot$x$data[[i]]$y <- frame_dt$count
          plot$x$data[[i]]$width <- 1
          plot$x$data[[i]]$marker$color <- frame_colors[[frame_name]]
          line_index <- i + length(frame_levels)
          plot$x$data[[line_index]]$x <- frame_dt$position
          plot$x$data[[line_index]]$y <- frame_dt$count
          plot$x$data[[line_index]]$line$color <- frame_colors[[frame_name]]
        }
      }
      keep_length <- if (is_gl_template && !has_switch_lines) {
        length(frame_levels)
      } else {
        length(frame_levels) * 2
      }
      if (length(plot$x$data) > keep_length) {
        plot$x$data <- plot$x$data[seq_len(keep_length)]
      }
    } else {
      for (frame_name in frame_levels) {
        plot <- plotly::add_trace(
          plot,
          x = c(NA_real_),
          y = c(NA_real_),
          text = NA_character_,
          type = "scattergl",
          mode = "lines",
          name = frame_name,
          legendgroup = frame_name,
          line = list(color = frame_colors[[frame_name]], width = 6, simplify = FALSE),
          hovertemplate = "%{text}<extra></extra>",
          meta = list(
            rc_columns_switch = "gl_subset",
            rc_columns_group = columns_trace_group,
            rc_columns_frame = frame_name
          ),
          showlegend = TRUE,
          visible = "legendonly",
          inherit = FALSE
        )
      }
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
          line = list(color = frame_colors[[frame_name]], width = 0.8, simplify = FALSE),
          hovertemplate = track_hover_template(include_frame = TRUE),
          meta = list(
            rc_columns_switch = "line",
            rc_columns_group = columns_trace_group,
            rc_columns_frame = frame_name
          ),
          showlegend = FALSE,
          visible = TRUE
        )
      }
    }
  } else if (type == "area") {
    frame_fill_colors <- track_frame_colors(frame_levels, frame_theme, with_alpha = TRUE)
    if (inherits(templates$cov_panel_area_plotly, "plotly")) {
      plot <- templates$cov_panel_area_plotly
      for (i in seq_along(frame_levels)) {
        frame_name <- frame_levels[[i]]
        frame_dt <- profile[as.character(frame) == frame_name]
        plot$x$data[[i]]$x <- frame_dt$position
        plot$x$data[[i]]$y <- frame_dt$count
        plot$x$data[[i]]$name <- frame_name
        plot$x$data[[i]]$legendgroup <- frame_name
        plot$x$data[[i]]$fillcolor <- grDevices::adjustcolor(
          frame_fill_colors[[frame_name]],
          alpha.f = track_area_fill_alpha()
        )
      }
      if (length(plot$x$data) > length(frame_levels)) {
        plot$x$data <- plot$x$data[seq_along(frame_levels)]
      }
    } else {
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
          line = list(color = "black", width = 0.4, simplify = FALSE),
          fill = "tozeroy",
          fillcolor = grDevices::adjustcolor(
            frame_fill_colors[[frame_name]],
            alpha.f = track_area_fill_alpha()
          ),
          hovertemplate = track_hover_template(include_frame = TRUE),
          showlegend = TRUE
        )
      }
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
        line = list(color = "black", width = 0.4, simplify = FALSE),
        fillcolor = grDevices::adjustcolor(frame_colors[[frame_name]], alpha.f = 0.8),
        hovertemplate = track_hover_template(include_frame = TRUE),
        showlegend = TRUE
      )
    }
    guide_y_max <- max(profile[, .(count = sum(count)), by = position]$count, na.rm = TRUE)
  }

  plot <- track_guides(plot, x_range, guide_y_max, lines, add_zero = TRUE, line_size = line_size)
  if (type == "columns") {
    plot <- ensure_columns_gl_hover(plot)
  }
  plot
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
