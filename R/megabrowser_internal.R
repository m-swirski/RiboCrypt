compute_collection_table_shiny <- function(mainPlotControls,
                                           path = mainPlotControls()$table_path,
                                           lib_sizes = mainPlotControls()$lib_sizes,
                                           df = mainPlotControls()$dff,
                                           metadata_field = mainPlotControls()$metadata_field,
                                           normalization = mainPlotControls()$normalization,
                                           kmer = mainPlotControls()$kmer,
                                           min_count = mainPlotControls()$min_count,
                                           subset = mainPlotControls()$subset,
                                           group_on_tx_tpm = mainPlotControls()$group_on_tx_tpm,
                                           split_by_frame = mainPlotControls()$frame,
                                           ratio_interval = mainPlotControls()$ratio_interval,
                                           enrichment_term = mainPlotControls()$enrichment_term,
                                           clusters = mainPlotControls()$clusters,
                                           metadata) {
  if (is.null(metadata)) stop("Metadata not defined, no metabrowser allowed for now!")
  time_before <- Sys.time()
  cat("Starting loading + Profile + plot calc\n")
  dtable <- compute_collection_table(path, lib_sizes, df, metadata_field,
                                     normalization, kmer, metadata, min_count,
                                     as_list = TRUE, subset = subset,
                                     group_on_tx_tpm = group_on_tx_tpm,
                                     split_by_frame = split_by_frame,
                                     ratio_interval = ratio_interval,
                                     enrichment_term = enrichment_term,
                                     clusters = clusters)
  timer_done_nice_print("Done: lib loading + Coverage calc: ", time_before)
  return(dtable)
}

mb_controller_shiny <- function(input, df, gene_name_list, cds, tx) {
  click_plot_browser_allsamp_controller(input, df, gene_name_list, cds, tx)
}

mb_plot_object_shiny <- function(table_obj, input, templates = NULL) {
  get_meta_browser_plot(table_obj, isolate(input$heatmap_color),
                        isolate(input$color_mult),
                        isolate(input$plotType),
                        template = templates$megabrowser_heatmap_plotly)
}

mb_mid_plot_shiny <- function(plot_object, plotType) {
  req(plotType == "plotly")
  plot_object$x$source <- "mb_mid"
  plotly::event_register(plot_object, "plotly_relayout")
}

mb_mid_image_shiny <- function(plot_object, session, ns, rel_heights = c(0.15, 0.75, 0.10)) {
  w <- session$clientData[[paste0("output_", ns("mb_subplot_static"), "_width")]]
  h <- session$clientData[[paste0("output_", ns("mb_subplot_static"), "_height")]]
  if (is.null(w) || is.na(w)) w <- 900
  if (is.null(h) || is.na(h)) h <- 700
  rel_heights <- rel_heights / sum(rel_heights)
  mid_h <- h * rel_heights[2]
  mid_w <- w
  plotly_image_from_plot(function() {
    if (is(plot_object, "Heatmap")) {
      pad_px <- c(t = 0, r = 0, b = 0, l = 0)
      pad <- grid::unit(pad_px, "mm")
      ComplexHeatmap::draw(plot_object, padding = pad)
    } else {
      print(plot_object)
    }
  }, width = mid_w, height = mid_h, res = 96)
}

mb_static_subplot_shiny <- function(top_plot, mid_image, bottom_plot, rel_heights = c(0.15, 0.75, 0.10)) {
  rel_heights <- rel_heights / sum(rel_heights)
  p <- suppressWarnings(plotly::subplot(
    top_plot,
    mid_image,
    bottom_plot,
    nrows = 3,
    shareX = TRUE,
    heights = rel_heights,
    margin = 0.01,
    titleX = FALSE,
    titleY = FALSE
  ))
  p %>%
    plotly::layout(
      dragmode = FALSE,
      xaxis = list(fixedrange = TRUE),
      xaxis2 = list(fixedrange = TRUE),
      xaxis3 = list(fixedrange = TRUE),
      yaxis = list(fixedrange = TRUE),
      yaxis2 = list(fixedrange = TRUE),
      yaxis3 = list(fixedrange = TRUE)
    ) %>%
    plotly::config(
      scrollZoom = FALSE,
      modeBarButtonsToRemove = c(
        "zoom2d", "pan2d", "select2d", "lasso2d",
        "zoomIn2d", "zoomOut2d", "autoScale2d", "resetScale2d"
      )
    )
}

mb_read_axis_event <- function(ed, axes) {
  if (!is.list(ed) && !is.environment(ed)) {
    return(list(r0 = NULL, r1 = NULL, r = NULL, auto = FALSE))
  }
  nms <- names(ed)
  if (is.null(nms) || length(nms) == 0) {
    return(list(r0 = NULL, r1 = NULL, r = NULL, auto = FALSE))
  }

  for (ax in axes) {
    k0 <- paste0(ax, ".range[0]")
    k1 <- paste0(ax, ".range[1]")
    kr <- paste0(ax, ".range")
    ka <- paste0(ax, ".autorange")
    if (!is.null(ed[[k0]]) && !is.null(ed[[k1]])) {
      return(list(r0 = ed[[k0]], r1 = ed[[k1]], r = NULL, auto = FALSE))
    }
    if (!is.null(ed[[kr]]) && length(ed[[kr]]) == 2) {
      return(list(r0 = NULL, r1 = NULL, r = ed[[kr]], auto = FALSE))
    }
    if (isTRUE(ed[[ka]])) {
      return(list(r0 = NULL, r1 = NULL, r = NULL, auto = TRUE))
    }
  }
  list(r0 = NULL, r1 = NULL, r = NULL, auto = FALSE)
}

sync_megabrowser_x_shiny <- function(ed, session, sync_tracks = TRUE, sync_sidebar = TRUE,
                                     y_max = NULL, y_reversed = TRUE,
                                     x_reset_range = NULL,
                                     x_axes = c("xaxis", "xaxis2", "xaxis3"),
                                     y_axes = c("yaxis", "yaxis2", "yaxis3")) {
  if (is.null(ed)) return(invisible(NULL))

  xed <- mb_read_axis_event(ed, x_axes)
  yed <- mb_read_axis_event(ed, y_axes)

  if (!is.null(xed$r0) && !is.null(xed$r1)) {
    relayout <- list("xaxis.range" = c(xed$r0, xed$r1), "xaxis.autorange" = FALSE)
  } else if (!is.null(xed$r) && length(xed$r) == 2) {
    relayout <- list("xaxis.range" = c(xed$r[[1]], xed$r[[2]]), "xaxis.autorange" = FALSE)
  } else if (xed$auto) {
    if (!is.null(x_reset_range) && length(x_reset_range) == 2) {
      relayout <- list("xaxis.range" = c(x_reset_range[[1]], x_reset_range[[2]]),
                       "xaxis.autorange" = FALSE)
    } else {
      relayout <- list("xaxis.autorange" = TRUE)
    }
  } else {
    relayout <- list()
  }

  if (sync_tracks) {
    plotly::plotlyProxy("mb_top_summary", session) %>%
      plotly::plotlyProxyInvoke("relayout", relayout)
    plotly::plotlyProxy("mb_bottom_gene", session) %>%
      plotly::plotlyProxyInvoke("relayout", relayout)
  }

  if (!is.null(yed$r0) && !is.null(yed$r1)) {
    y0_raw <- yed$r0
    y1_raw <- yed$r1
    if (y_reversed && !is.null(y_max)) {
      y0_raw <- (y_max + 1) - y0_raw
      y1_raw <- (y_max + 1) - y1_raw
    }
    y0 <- min(y0_raw, y1_raw)
    y1 <- max(y0_raw, y1_raw)
    y0c <- round(y0)
    y1c <- round(y1)
    if (y0c > y1c) {
      tmp <- y0c
      y0c <- y1c
      y1c <- tmp
    }
    y0e <- y0c - 0.5
    y1e <- y1c + 0.5
    yrelayout <- list(
      yaxis = list(
        range = c(y1e, y0e), autorange = FALSE,
        showticklabels = FALSE, ticks = "", showgrid = FALSE, zeroline = FALSE, title = NULL
      )
    )
  } else if (!is.null(yed$r) && length(yed$r) == 2) {
    y0_raw <- yed$r[[1]]
    y1_raw <- yed$r[[2]]
    if (y_reversed && !is.null(y_max)) {
      y0_raw <- (y_max + 1) - y0_raw
      y1_raw <- (y_max + 1) - y1_raw
    }
    y0 <- min(y0_raw, y1_raw)
    y1 <- max(y0_raw, y1_raw)
    y0c <- round(y0)
    y1c <- round(y1)
    if (y0c > y1c) {
      tmp <- y0c
      y0c <- y1c
      y1c <- tmp
    }
    y0e <- y0c - 0.5
    y1e <- y1c + 0.5
    yrelayout <- list(
      yaxis = list(
        range = c(y1e, y0e), autorange = FALSE,
        showticklabels = FALSE, ticks = "", showgrid = FALSE, zeroline = FALSE, title = NULL
      )
    )
  } else if (yed$auto) {
    yrelayout <- list(
      yaxis = list(
        autorange = "reversed",
        showticklabels = FALSE, ticks = "", showgrid = FALSE, zeroline = FALSE, title = NULL
      )
    )
  } else {
    yrelayout <- NULL
  }

  if (sync_sidebar && !is.null(yrelayout)) {
    plotly::plotlyProxy("d", session) %>%
      plotly::plotlyProxyInvoke("relayout", yrelayout)
  }
  invisible(NULL)
}

megabrowser_full_x_range <- function(display_range = NULL, table = NULL) {
  if (!is.null(display_range)) {
    max_pos <- suppressWarnings(as.numeric(widthPerGroup(display_range, FALSE)))
    if (length(max_pos) > 1) max_pos <- max_pos[[1]]
    if (is.finite(max_pos) && max_pos >= 1) {
      return(c(1, max_pos))
    }
  }
  summary_cov <- attr(table, "summary_cov")
  if (is.data.frame(summary_cov) &&
      nrow(summary_cov) > 0 &&
      "position" %in% names(summary_cov)) {
    max_pos <- suppressWarnings(max(summary_cov$position, na.rm = TRUE))
    if (is.finite(max_pos) && max_pos >= 1) {
      return(c(1, max_pos))
    }
  }
  c(1, nrow(table))
}

get_meta_browser_plot <- function(table, color_theme, color_mult = 3,
                                  plotType = "plotly", template = NULL) {
  time_before <- Sys.time()
  cat("Creating metabrowser heatmap\n")

  colors <- if (color_theme == "default (White-Blue)") {
    c("white", "lightblue", rep("blue", 4 + color_mult), "navy", "black")
  } else if (color_theme == "Matrix (black,green,red)") {
    c("#000000", "#2CFA1F", "yellow2", rep("#FF2400", color_mult))
  } else stop("Invalid color theme!")
  mat <- t(table)

  km <- attr(table, "km")
  stopifnot(!is.null(km))
  row_clusters <- attr(table, "row_order_list")
  stopifnot(!is.null(row_clusters))

  if (plotType == "plotly") {
    mat <- mat[rev(unlist(row_clusters, use.names = FALSE)),]
    plot_mat <- mat[rev(seq_len(nrow(mat))), , drop = FALSE]
    ratio <- attr(table, "ratio")
    if (is.null(ratio) || !is.numeric(ratio) || length(ratio) != 1 || !is.finite(ratio) || ratio < 1) {
      ratio <- 1
    }
    ratio <- as.integer(ratio)
    x_positions <- ((seq_len(ncol(mat)) - 1L) * ratio) + 1L
    full_x_range <- megabrowser_full_x_range(table = table)
    full_y_range <- c(0.5, nrow(mat) + 0.5)
    colorscale <- megabrowser_heatmap_colorscale(colors)

    margins <- margin_megabrowser()
    margins$t <- margins$b  <- 8
    if (inherits(template, "plotly")) {
      plot <- template
      plot$x$data[[1]]$x <- x_positions
      plot$x$data[[1]]$y <- NULL
      plot$x$data[[1]]$z <- plot_mat
      plot$x$data[[1]]$colorscale <- colorscale
      plot$x$data[[1]]$zmin <- NULL
      plot$x$data[[1]]$zmax <- NULL
      plot$x$data[[1]]$zauto <- NULL
      plot$x$data[[1]]$autocolorscale <- NULL
      if (length(plot$x$attrs) >= 1) {
        plot$x$attrs[[1]]$x <- x_positions
        plot$x$attrs[[1]]$y <- NULL
        plot$x$attrs[[1]]$z <- plot_mat
        plot$x$attrs[[1]]$colorscale <- colorscale
      }
      plot <- plot %>% plotly::layout(
        margin = margins,
        xaxis = list(
          range = full_x_range,
          autorange = FALSE
        ),
        yaxis = list(
          range = full_y_range,
          autorange = FALSE
        ),
        dragmode = "zoom"
      ) %>%
        plotly::config(doubleClick = FALSE)
    } else {
      plot <- plotly::plot_ly(
        x = x_positions,
        z = plot_mat,
        colors = colors,
        showscale = FALSE,
        type = "heatmapgl"
      ) %>% plotly::layout(
        margin = margins,
        xaxis = list(
          range = full_x_range,
          autorange = FALSE
        ),
        yaxis = list(
          range = full_y_range,
          autorange = FALSE
        ),
        dragmode = "zoom"
      ) %>%
        plotly::config(doubleClick = FALSE)
    }

  } else {
    mat <- mat[rev(unlist(row_clusters, use.names = FALSE)),]
    cluster <- km$cluster[rev(unlist(row_clusters, use.names = FALSE))]
    cluster <- factor(cluster, levels = rev(sort(unique(cluster))))

    plot <- ComplexHeatmap::Heatmap(mat, show_row_dend = FALSE,
                            cluster_columns = FALSE,
                            cluster_rows = FALSE,
                            use_raster = TRUE,  raster_quality = 5,
                            split = cluster, gap = unit(0.2, "mm"),
                            col =  colors, show_row_names = FALSE,
                            show_heatmap_legend = FALSE, row_title = NULL)
  }
  timer_done_nice_print("-- metabrowser heatmap done: ", time_before)
  return(plot)
}

megabrowser_heatmap_colorscale <- function(colors) {
  rgb <- grDevices::col2rgb(colors)
  rgba <- grDevices::rgb(rgb[1, ], rgb[2, ], rgb[3, ], maxColorValue = 255)
  cbind(
    seq(0, 1, length.out = length(rgba)),
    rgba
  )
}

megabrowserHeatmapPlotlyTemplate <- function() {
  plotly::plotly_build(
    plotly::plot_ly(
      x = 1,
      y = 1,
      z = matrix(0, nrow = 1),
      type = "heatmapgl",
      showscale = FALSE,
      hoverinfo = "skip"
    ) %>%
      plotly::layout(
        margin = margin_megabrowser(),
        xaxis = list(autorange = FALSE),
        yaxis = list(autorange = FALSE),
        dragmode = "zoom"
      )
  )
}

get_meta_browser_plot_full_shiny <- function(table, plot_object, controller, theme_template) {
  get_meta_browser_plot_full(table()$table,
                             plot_object(), controller()$id,
                             controller()$dff, controller()$summary_track,
                             controller()$display_annot,
                             region_type = controller()$region_type,
                             controller()$plotType, controller()$tx_annotation,
                             controller()$display_region,
                             controller()$annotation,
                             controller()$viewMode,
                             controller()$collapsed_introns_width,
                             theme_template
  )
}

renderMegabrowser <- function(plotType, ns,
                              height = "700px", width = "100%",
                              rel_heights = c(0.15, 0.75, 0.10)) {
  time_before <- Sys.time()
  # split container height by rel_heights
  rel_heights <- rel_heights / sum(rel_heights)
  h1 <- sprintf("calc(%s * %.6f)", height, rel_heights[1])
  h2 <- sprintf("calc(%s * %.6f)", height, rel_heights[2])
  h3 <- sprintf("calc(%s * %.6f)", height, rel_heights[3])

  if (plotType == "ggplot2") {
    tagList <- plotly::plotlyOutput(ns("mb_subplot_static"), height = height, width = width) %>%
      shinycssloaders::withSpinner(color = "#0dc5c1")
  } else {
    top <- plotly::plotlyOutput(ns("mb_top_summary"), height = h1, width = width)
    middle <- plotly::plotlyOutput(ns("myPlotlyPlot"), height = h2, width = width)
    middle <- middle %>% shinycssloaders::withSpinner(color = "#0dc5c1")
    bottom <- plotly::plotlyOutput(ns("mb_bottom_gene"), height = h3, width = width)
    tagList <- tagList(top, middle, bottom)
  }
  timer_done_nice_print("-- Mega browser 3 row plot done: ", time_before)
  return(tagList)
}

plotly_image_from_plot <- function(plot_fn, width, height, res = 96,
                                   bg = "transparent") {
  if (!requireNamespace("ragg", quietly = TRUE)) {
    stop("Package 'ragg' is required for plotly image rendering.")
  }
  if (!requireNamespace("base64enc", quietly = TRUE)) {
    stop("Package 'base64enc' is required for plotly image rendering.")
  }
  tmp <- tempfile(fileext = ".png")
  ragg::agg_png(filename = tmp, width = width, height = height,
                units = "px", res = res, background = bg)
  on.exit({
    if (file.exists(tmp)) unlink(tmp)
  }, add = TRUE)
  plot_fn()
  grDevices::dev.off()
  data_uri <- base64enc::dataURI(file = tmp, mime = "image/png")
  plotly::plot_ly(
    x = c(0, 1),
    y = c(0, 1),
    type = "scatter",
    mode = "markers",
    marker = list(opacity = 0),
    hoverinfo = "skip",
    showlegend = FALSE
  ) %>%
    plotly::layout(
      images = list(list(
        source = data_uri,
        xref = "paper", yref = "paper",
        x = 0, y = 1, sizex = 1, sizey = 1,
        sizing = "stretch", layer = "below"
      )),
      xaxis = list(visible = FALSE, range = c(0, 1), fixedrange = TRUE),
      yaxis = list(visible = FALSE, range = c(0, 1), fixedrange = TRUE),
      margin = list(l = 0, r = 0, t = 0, b = 0),
      paper_bgcolor = "rgba(0,0,0,0)",
      plot_bgcolor  = "rgba(0,0,0,0)"
    )
}

get_megabrowser_annotation_plot_shiny <- function(controller, templates = NULL) {
  p <- get_megabrowser_annotation_plot(controller()$id,
                             controller()$dff, controller()$summary_track,
                             controller()$display_annot,
                             region_type = controller()$region_type,
                             controller()$plotType, controller()$tx_annotation,
                             controller()$display_region,
                             controller()$annotation,
                             controller()$customRegions,
                             controller()$viewMode,
                             controller()$collapsed_introns_width,
                             templates = templates
  )
  p$x$source <- "mb_bottom"
  p
}

#' Full plot for allsamples browser
#'
#' @param m data.table of coverage per sample (wide format)
#' @param heatmap ComplexHeatmap object of plot from 'm'
#' @param id id of transcript
#' @param df ORFik experiment
#' @param summary logical, default TRUE (add top plot)
#' @param annotation logical, default TRUE (add bottom annotation track)
#' @param region_type character, "what is the coverage region?" Usually full mrna:
#' "mrna" or "leader+cds".
#' @param plotType = "plotly",
#' @param tx_annotation a GRangesList of tx annotation
#' @param display_region a GRangesList of display region
#' @param cds_annotation a GRangesList of cds annotation
#' @param viewMode character, "tx" or "genomic"
#' @param collapse_intron_flank integer, if TRUE and viewMode genomic, collapse
#' introns to this max size.
#' @param rel_heights numeric < 1, default: c(0.2, 0.75, 0.05).
#' Relative heights, sum to 1 and must be length 3.
#' @return a cowplot grub
#' @importFrom cowplot plot_grid
get_megabrowser_annotation_plot <- function(id, df,
                                       summary = TRUE, annotation = TRUE,
                                       region_type,
                                       plotType = "plotly",
                                       tx_annotation, display_region,
                                       cds_annotation, custom_regions = NULL,
                                       viewMode,
                                       collapse_intron_flank,
                                       templates = NULL) {
  time_before <- Sys.time()
  res <- annotation_track_allsamples(df, id, display_region, cds_annotation,
                                     tx_annotation, custom_regions,
                                     viewMode, collapse_intron_flank,
                                     templates = templates)

  timer_done_nice_print("-- Mega browser annotation plot done: ", time_before)
  return(res)
}

add_alpha <- function(cols, alpha = 0.7) {
  rgb <- grDevices::col2rgb(cols)
  sprintf("rgba(%d,%d,%d,%.3f)",
          rgb[1, ], rgb[2, ], rgb[3, ], alpha)
}

profile_plotly_gl <- function(dt, frame_colors = frame_color_themes("R"),
                              bar_px = 6, alpha = 0.7,
                              template = covPanelColumnsGLPlotlyTemplate(bar_px = bar_px)) {
  # dt needs: position, count, frame
  stopifnot(is.data.table(dt) || is.data.frame(dt))
  dt <- as.data.table(dt)
  x_range <- range(dt$position, na.rm = TRUE)


  # ensure deterministic order
  dt[, frame := as.character(frame)]
  frame_colors <- add_alpha(frame_colors, alpha)

  p <- template
  frame_levels <- sort(unique(dt$frame))
  for (i in seq_along(frame_levels)) {
    fr <- frame_levels[[i]]
    d <- dt[frame == fr]
    if (nrow(d) == 0L) next

    tt <- paste0(
      "frame: ", fr,
      "<br>pos: ", d$position,
      "<br>count: ", d$count
    )
    seg <- track_vertical_segments(d$position, d$count, tt)

    col <- frame_colors[as.integer(fr) + 1L]
    if (is.null(col)) col <- "grey50"

    p$x$data[[i]]$x <- seg$x
    p$x$data[[i]]$y <- seg$y
    p$x$data[[i]]$text <- seg$text
    p$x$data[[i]]$name <- fr
    p$x$data[[i]]$legendgroup <- fr
    p$x$data[[i]]$line$color <- col
  }
  if (length(p$x$data) > length(frame_levels)) {
    p$x$data <- p$x$data[seq_along(frame_levels)]
  }

  p %>%
    layout(
      showlegend = FALSE,        # remove legend globally
      margin = margin_megabrowser(),
      paper_bgcolor = "rgba(0,0,0,0)",
      plot_bgcolor  = "rgba(0,0,0,0)",
      xaxis = list(
        range = x_range,
        autorange = FALSE,
        showticklabels = FALSE,
        ticks = "",
        showgrid = FALSE,
        zeroline = FALSE,
        title = NULL,
        fixedrange = FALSE
      ),
      yaxis = list(
        showticklabels = FALSE,
        ticks = "",
        showgrid = FALSE,
        zeroline = FALSE,
        title = NULL,
        fixedrange = TRUE
      )
    )
}

summary_track_allsamples <- function(mat, template = NULL) {
  time_before <- Sys.time()
  if (is.null(template)) {
    template <- covPanelColumnsPlotlyTemplate()
  }
  is_empty_summary <- is.null(mat) ||
    nrow(mat) == 0L ||
    !all(c("position", "count") %in% names(mat)) ||
    !any(is.finite(mat$position)) ||
    !any(is.finite(mat$count))

  if (is_empty_summary) {
    summary_plot <- template %>%
      plotly::layout(
        showlegend = FALSE,
        hovermode = "closest",
        margin = margin_megabrowser(),
        paper_bgcolor = "rgba(0,0,0,0)",
        plot_bgcolor = "rgba(0,0,0,0)",
        xaxis = list(
          range = c(0, 1),
          autorange = FALSE,
          showticklabels = FALSE,
          ticks = "",
          showgrid = FALSE,
          zeroline = FALSE,
          title = NULL,
          fixedrange = FALSE
        ),
        yaxis = list(
          showticklabels = FALSE,
          ticks = "",
          showgrid = FALSE,
          zeroline = FALSE,
          title = NULL,
          fixedrange = TRUE
        )
      ) %>%
      plotly::config(doubleClick = FALSE)
    summary_plot$x$layout$yaxis$title <- NULL
    summary_plot$x$source <- "mb_top"
    timer_done_nice_print("-- Mega browser summary track done: ", time_before)
    return(summary_plot)
  }

  x_range <- range(mat$position, na.rm = TRUE)
  summary_plot <- createSinglePlot(
    profile = mat,
    withFrames = TRUE,
    frame_colors = "R",
    colors = NULL,
    ylabels = "summary",
    ylabels_full_name = "summary",
    lines = numeric(),
    type = "columns",
    lib_index = 1,
    total_libs = 1,
    flip_ylabel = FALSE,
    templates = list(cov_panel_columns_plotly = template)
  ) %>%
    plotly::layout(
      showlegend = FALSE,
      hovermode = "closest",
      margin = margin_megabrowser(),
      paper_bgcolor = "rgba(0,0,0,0)",
      plot_bgcolor  = "rgba(0,0,0,0)",
      xaxis = list(
        range = x_range,
        autorange = FALSE,
        showticklabels = FALSE,
        ticks = "",
        showgrid = FALSE,
        zeroline = FALSE,
        title = NULL,
        fixedrange = FALSE
      ),
      yaxis = list(
        showticklabels = FALSE,
        ticks = "",
        showgrid = FALSE,
        zeroline = FALSE,
        title = NULL,
        fixedrange = TRUE
      )
    ) %>%
    addColumnsZoomSwitch() %>%
    lineDeSimplify() %>%
    plotly::config(doubleClick = FALSE)
  summary_plot$x$layout$yaxis$title <- NULL
  summary_plot$x$source <- "mb_top"
  timer_done_nice_print("-- Mega browser summary track done: ", time_before)
  return(summary_plot)
}

annotation_track_allsamples <- function(df, id, display_range, annotation,
                                        tx_annotation, custom_regions = NULL,
                                        viewMode, collapse_intron_flank,
                                        templates = NULL) {
  annotation_list <- annotation_controller(
    df = df,
    display_range = display_range,
    annotation = annotation,
    leader_extension = 0,
    trailer_extension = 0,
    viewMode = viewMode
  )
  gene_model_panel <- createGeneModelPanel(annotation_list$display_range, annotation_list$annotation,
                                           tx_annotation = tx_annotation,
                                           custom_regions = custom_regions,
                                           viewMode = viewMode, collapse_intron_flank)
  return(geneModelPanelPlotly(
    gene_model_panel[[1]],
    template = templates$gene_model_panel_plotly
  ))
}

get_lib_sizes_file <- function(dff) {
  lib_sizes <- file.path(QCfolder(dff), "totalCounts_mrna.rds")
  if (!file.exists(lib_sizes))
    stop("Count table library size files are not created, missing file totalCounts_mrna.rds",
         " see vignette for more information on how to make these.")
  return(lib_sizes)
}

validate_enrichment_term <- function(enrichment_term, clusters, ratio_interval, other_gene, metadata_field) {
  valid_enrichment_clusterings <- c(clusters != 1, isTruthy(ratio_interval), isTruthy(other_gene))
  enrichment_test_types <- c("Clusters", "Ratio bins", "Other gene tpm bins")[valid_enrichment_clusterings]

  valid_enrichment_terms <- c(metadata_field, enrichment_test_types)
  shiny::validate(
    shiny::need(
      enrichment_term %in% valid_enrichment_terms,
      paste0(
        "Enrichment term is not valid. Valid options: ",
        paste(valid_enrichment_terms, collapse = ", ")
      )
    )
  )
  invisible(NULL)
}

#' Given ratio interval user input, format correctly
#' @param ratio_interval character or NULL, if defined character, will format it
#' @return a numeric interval, if input !isTruthy, returns NULL
#' @noRd
get_ratio_interval <- function(ratio_interval) {
  if (isTruthy(ratio_interval)) {
    split_ratio <- unlist(strsplit(ratio_interval, ";"))

    temp_interval <- strsplit(split_ratio, ":|-")
    single_input <- lengths(temp_interval) == 1
    if (any(lengths(temp_interval) > 2)) stop("Malformed sort in interval input!")
    if(any(single_input)) {
      temp_interval[which(single_input)] <- lapply(temp_interval[which(single_input)], function(x) rep(x, 2))
    }
    temp_interval <- as.numeric(unlist(temp_interval))

    if (anyNA(temp_interval)) stop("Malformed sort in interval input!")
    ratio_interval <- temp_interval
  } else ratio_interval <- NULL

  return(ratio_interval)
}
