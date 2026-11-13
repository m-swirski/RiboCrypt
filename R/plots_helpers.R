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
    # WebGL column traces, or column traces drawn as SVG on phones (marked gl_subset); unmarked SVG lines
    # (guides, smoothed lines) keep their own hover
    is_columns_gl_trace <- identical(trace$mode, "lines") &&
      (identical(trace$type, "scattergl") || identical(trace$meta$rc_columns_switch, "gl_subset"))
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

# Plotly keeps every original point for hover and data export; this only
# simplifies SVG path drawing when a trace greatly exceeds the display width.
browser_simplify_coverage <- function(n) isTRUE(getOption("RiboCrypt.browser_simplify", TRUE)) && n > 5000L

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
                                        templates = NULL, is_cellphone = FALSE) {
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
      plot$x$data[[1]]$line$simplify <- browser_simplify_coverage(nrow(profile))
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
        line = list(color = "black", width = 0.5,
                    simplify = type %in% c("lines", "area") && browser_simplify_coverage(nrow(profile))),
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
        plot$x$data[[i]]$line$simplify <- browser_simplify_coverage(nrow(frame_dt))
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
          line = list(color = frame_colors[[frame_name]], width = 0.8,
                      simplify = browser_simplify_coverage(nrow(frame_dt))),
          hovertemplate = track_hover_template(include_frame = TRUE),
          showlegend = TRUE
        )
      }
    }
  } else if (type == "columns") {
    if (inherits(templates$cov_panel_columns_plotly, "plotly")) {
      plot <- templates$cov_panel_columns_plotly
      is_gl_template <- identical(plot$x$data[[1]]$type, "scattergl") ||
        identical(plot$x$data[[1]]$meta$rc_columns_switch, "gl_subset")
      if (is_gl_template) {
        for (i in seq_len(min(length(frame_levels), length(plot$x$data))))
          plot$x$data[[i]]$type <- if (isTRUE(is_cellphone)) "scatter" else "scattergl"
      }
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
          type = if (isTRUE(is_cellphone)) "scatter" else "scattergl",
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
        plot$x$data[[i]]$line$simplify <- browser_simplify_coverage(nrow(frame_dt))
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
          line = list(color = "black", width = 0.4,
                      simplify = browser_simplify_coverage(nrow(frame_dt))),
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
    plot <- ensure_columns_gl_hover(plotly::plotly_build(plot))
  }
  plot
}

#' Normalize the number of coverage tracks used for axis styling.
#' @noRd
browser_library_count <- function(total_libs) {
  if (!length(total_libs)) return(1L)
  count <- suppressWarnings(as.integer(total_libs[[1]]))
  if (is.na(count) || count < 1L) 1L else count
}

#' Number of non-zero y-axis labels that fit in a browser coverage track.
#' @noRd
browser_y_tick_count <- function(total_libs) {
  total_libs <- browser_library_count(total_libs)
  if (total_libs <= 2L) return(3L)
  if (total_libs == 3L) return(2L)
  1L
}

#' Font size for browser coverage y-axis labels.
#' @noRd
browser_y_tick_font_size <- function(total_libs) {
  total_libs <- browser_library_count(total_libs)

  if (total_libs <= 1L) return(14L)
  if (total_libs <= 3L) return(13L)
  if (total_libs <= 8L) return(11L)
  if (total_libs <= 12L) return(10L)
  if (total_libs <= 20L) return(9L)
  if (total_libs <= 35L) return(8L)
  if (total_libs <= 50L) return(7L)
  6L
}

#' Font size for browser coverage track names.
#' @noRd
browser_y_title_font_size <- function(total_libs) {
  total_libs <- browser_library_count(total_libs)

  if (total_libs <= 1L) return(16L)
  if (total_libs <= 2L) return(14L)
  if (total_libs <= 4L) return(12L)
  if (total_libs <= 8L) return(10L)
  if (total_libs <= 20L) return(9L)
  if (total_libs <= 35L) return(8L)
  if (total_libs <= 50L) return(7L)
  6L
}

#' Find a readable step close to the requested numeric interval.
#' @noRd
browser_nice_y_step <- function(target, round_up = FALSE) {
  if (length(target) != 1L || !is.finite(target) || target <= 0) return(1)

  exponent <- floor(log10(target))
  multipliers <- c(1, 1.25, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10)
  candidates <- as.vector(outer(multipliers, 10 ^ ((exponent - 1L):(exponent + 1L))))
  candidates <- unique(candidates[is.finite(candidates) & candidates > 0])
  if (isTRUE(round_up)) candidates <- candidates[candidates >= target]
  distance <- abs(candidates - target)

  # On an exact tie, prefer the value above the target. A slightly generous
  # scale reads better than a tick that substantially understates the peak.
  candidates[order(distance, candidates < target, candidates)][[1]]
}

#' Compact, deterministic labels for browser coverage ticks.
#' @noRd
browser_format_y_tick <- function(value) {
  if (!length(value)) return("")
  value <- as.numeric(value[[1]])
  if (!is.finite(value)) return("")

  magnitude <- abs(value)
  scale <- 1
  suffix <- ""
  if (magnitude >= 1e12) {
    scale <- 1e12
    suffix <- "T"
  } else if (magnitude >= 1e9) {
    scale <- 1e9
    suffix <- "B"
  } else if (magnitude >= 1e6) {
    scale <- 1e6
    suffix <- "M"
  } else if (magnitude >= 1e3) {
    scale <- 1e3
    suffix <- "k"
  }

  scaled <- signif(value / scale, 3)
  if (scale == 1 && abs(scaled) > 0 && abs(scaled) < 1e-3) {
    label <- format(scaled, scientific = TRUE, trim = TRUE, digits = 3)
    label <- sub("e([+-])0+", "e\\1", label)
  } else {
    label <- format(scaled, scientific = FALSE, trim = TRUE, digits = 3)
  }
  if (grepl(".", label, fixed = TRUE)) {
    label <- sub("0+$", "", label)
    label <- sub("[.]$", "", label)
  }
  paste0(label, suffix)
}

#' Convert proposed coverage ticks to positive integers where counts require it.
#' @noRd
browser_integer_y_ticks <- function(values, y_max, tick_count) {
  maximum_integer <- max(1, floor(y_max))
  desired_count <- as.integer(min(tick_count, maximum_integer))
  values <- sort(unique(pmin(
    maximum_integer,
    pmax(1, floor(values))
  )))

  # Flooring can merge small ticks (for example 0.8, 1.6 and 2.4). Fill the
  # resulting gaps when the count range contains enough distinct integers.
  if (length(values) < desired_count) {
    candidate <- 1
    while (length(values) < desired_count && candidate <= maximum_integer) {
      if (!candidate %in% values) values <- c(values, candidate)
      candidate <- candidate + 1
    }
    values <- sort(values)
  }
  values
}

#' Explicit non-zero coverage ticks and a range that always contains them.
#' @noRd
browser_y_ticks <- function(y_max, total_libs, integer_only = FALSE) {
  if (!length(y_max)) y_max <- 0
  y_max <- suppressWarnings(as.numeric(y_max[[1]]))
  if (!length(y_max) || !is.finite(y_max) || y_max <= 0) {
    return(list(values = 1, text = "1", range = c(0, 1.08)))
  }

  tick_count <- browser_y_tick_count(total_libs)
  target_fraction <- c(`1` = 0.82, `2` = 0.84, `3` = 0.88)[[as.character(tick_count)]]
  target_step <- y_max * target_fraction / tick_count
  step <- browser_nice_y_step(target_step)
  minimum_step <- y_max * 0.78 / tick_count
  if (step < minimum_step) {
    step <- browser_nice_y_step(minimum_step, round_up = TRUE)
  }
  values <- signif(seq_len(tick_count) * step, 12)
  if (isTRUE(integer_only)) {
    values <- browser_integer_y_ticks(values, y_max, tick_count)
  }
  range_top <- max(y_max, values) * 1.08

  list(
    values = values,
    text = vapply(values, browser_format_y_tick, character(1)),
    range = c(0, range_top)
  )
}

#' Whether the rendered coverage values can be labelled as integer counts.
#' @noRd
browser_track_has_integer_counts <- function(profile) {
  profile <- data.table::as.data.table(profile)
  if (!"count" %in% names(profile)) return(FALSE)

  count_values <- suppressWarnings(as.numeric(profile$count))
  count_values <- count_values[is.finite(count_values)]
  if (!length(count_values)) return(FALSE)

  tolerance <- sqrt(.Machine$double.eps) * pmax(1, abs(count_values))
  all(abs(count_values - round(count_values)) <= tolerance)
}

#' Maximum actually rendered by a browser coverage track.
#' @noRd
browser_track_y_max <- function(profile, type) {
  profile <- data.table::as.data.table(profile)
  if (!"count" %in% names(profile)) return(0)

  count_values <- suppressWarnings(as.numeric(profile$count))
  finite <- is.finite(count_values)
  if (!any(finite)) return(0)

  if (identical(type, "stacks") && "position" %in% names(profile)) {
    count_values <- as.numeric(rowsum(
      count_values[finite],
      group = profile$position[finite],
      reorder = FALSE
    ))
  } else {
    count_values <- count_values[finite]
  }

  max(c(0, count_values), na.rm = TRUE)
}

#' Coverage y-axis shared by the regular and observatory browsers.
#' @noRd
browser_coverage_yaxis <- function(y_max, total_libs, title = "",
                                   integer_only = FALSE) {
  ticks <- browser_y_ticks(y_max, total_libs, integer_only = integer_only)
  list(
    autorange = FALSE,
    range = ticks$range,
    # Keep the deterministic initial range without disabling Plotly's
    # interactive y-axis zoom.
    fixedrange = FALSE,
    zeroline = TRUE,
    zerolinecolor = "rgba(0,0,0,0.45)",
    zerolinewidth = 0.75,
    showgrid = TRUE,
    gridcolor = "rgba(0,0,0,0.10)",
    gridwidth = 0.5,
    tickmode = "array",
    tickvals = ticks$values,
    ticktext = ticks$text,
    ticks = "",
    showticklabels = TRUE,
    tickfont = list(
      size = browser_y_tick_font_size(total_libs),
      color = "#444444"
    ),
    title = list(
      text = title,
      font = list(size = browser_y_title_font_size(total_libs)),
      standoff = 8
    )
  )
}

singlePlot_add_theme <- function(profile_plot, ylabels, type,
                                 flip_ylabel = type == "heatmap", total_libs,
                                 ylabels_full_name = ylabels, as_plotly = TRUE,
                                 y_max = 0, integer_ticks = FALSE) {
  total_libs <- browser_library_count(total_libs)
  y_text_size <- browser_y_title_font_size(total_libs)
  annotation_list <- list()

  yaxis <- if (type == "heatmap") {
    list(
      autorange = TRUE,
      rangemode = "normal",
      zeroline = FALSE,
      showticklabels = FALSE,
      ticks = "",
      title = list(text = "")
    )
  } else {
    browser_coverage_yaxis(
      y_max, total_libs, ylabels,
      integer_only = integer_ticks
    )
  }

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
    layout_args$yaxis$title <- list(text = "")

    if (flip_ylabel || total_libs > 5) {
      dense_index_label <- total_libs > 5L
      # Every dense subplot uses the same pixel anchor. In particular, this is
      # deliberately independent of whether its tick reads "2" or "2.5".
      annotation_shift <- if (dense_index_label) -52 else -10
      annotation_characters <- if (dense_index_label) {
        nchar(as.character(total_libs))
      } else {
        nchar(as.character(ylabels))
      }
      annotation_width <- annotation_characters * y_text_size * 0.58
      layout_args$margin$l <- ceiling(min(160, max(45,
        abs(annotation_shift) + annotation_width + 6
      )))
      annotation_list <- list(list(
        text = ylabels,
        x = 0,
        y = 0.5,
        xref = "paper",
        yref = "y domain",
        xanchor = "right",
        yanchor = "middle",
        xshift = annotation_shift,
        showarrow = FALSE,
        font = list(size = y_text_size),
        hovertext = ylabels_full_name
      ))
    }
  }

  if (length(annotation_list) > 0) layout_args$annotations <- annotation_list

  do.call(plotly::layout, c(list(profile_plot), layout_args))
}
