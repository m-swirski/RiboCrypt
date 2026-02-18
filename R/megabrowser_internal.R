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
                                           metadata) {
  if (is.null(metadata)) stop("Metadata not defined, no metabrowser allowed for now!")
  time_before <- Sys.time()
  cat("Starting loading + Profile + plot calc\n")
  # browser()
  dtable <- compute_collection_table(path, lib_sizes, df, metadata_field,
                                     normalization, kmer, metadata, min_count,
                                     as_list = TRUE, subset = subset,
                                     group_on_tx_tpm = group_on_tx_tpm,
                                     split_by_frame = split_by_frame,
                                     ratio_interval = ratio_interval,
                                     enrichment_term = enrichment_term)
  timer_done_nice_print("Done: lib loading + Coverage calc: ", time_before)
  return(dtable)
}

get_meta_browser_plot <- function(table, color_theme, clusters = 1,
                                  color_mult = 3, plotType = "plotly") {
  time_before <- Sys.time()
  cat("Creating metabrowser heatmap\n")

  colors <- if (color_theme == "default (White-Blue)") {
    c("white", "lightblue", rep("blue", 4 + color_mult), "navy", "black")
  } else if (color_theme == "Matrix (black,green,red)") {
    c("#000000", "#2CFA1F", "yellow2", rep("#FF2400", color_mult))
  } else stop("Invalid color theme!")
  mat <- t(table)

  km <- kmeans(mat, centers = clusters)
  row_clusters <- seq_along(km$cluster)
  row_clusters <- split(row_clusters, km$cluster)


  if (plotType == "plotly") {
    mat <- mat[unlist(row_clusters, use.names = FALSE),]
    ratio <- attr(table, "ratio")
    x_range <- which(!duplicated((seq_len(nrow(mat)) - 1L) %/% ratio + 1L))
    plot <- plotly::plot_ly(
      x = ~x_range,
      z = mat,
      colors = colors,
      showscale = FALSE,
      type = "heatmapgl"
    ) %>% plotly::layout(margin = list(l = 30, r = 100, t = 10, b = 82)) %>%
      plotly::config(doubleClick = "reset")
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
  attr(plot, "row_order_list") <- row_clusters
  attr(plot, "clusters") <- length(row_clusters)

  timer_done_nice_print("-- metabrowser heatmap done: ", time_before)
  return(plot)
}

row_order <- function(x) {
  if (!is.null(attr(x, "row_order_list"))) {
    attr(x, "row_order_list")
  } else if (is(x, "Heatmap")) {
    ComplexHeatmap::row_order(x)
  } else stop("row_order must be attribute or x must be ComplexHeatmap!")
}

get_meta_browser_plot_full_shiny <- function(table, plot_object, controller, gg_theme) {
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
                             gg_theme
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

  top <- plotly::plotlyOutput(ns("mb_top_summary"), height = h1, width = width)

  middle <- if (plotType == "plotly") {
    plotly::plotlyOutput(ns("myPlotlyPlot"), height = h2, width = width)
  } else {
    plotOutput(ns("myGgplot"), height = h2, width = width)
  }

  middle <- middle %>% shinycssloaders::withSpinner(color = "#0dc5c1")

  bottom <- plotly::plotlyOutput(ns("mb_bottom_gene"), height = h3, width = width)
  tagList <- tagList(top, middle, bottom)
  timer_done_nice_print("-- Mega browser 3 row plot done: ", time_before)
  return(tagList)
}

get_megabrowser_annotation_plot_shiny <- function(controller, gg_theme) {
  get_megabrowser_annotation_plot(controller()$id,
                             controller()$dff, controller()$summary_track,
                             controller()$display_annot,
                             region_type = controller()$region_type,
                             controller()$plotType, controller()$tx_annotation,
                             controller()$display_region,
                             controller()$annotation,
                             controller()$viewMode,
                             controller()$collapsed_introns_width,
                             gg_theme
  )
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
                                       cds_annotation, viewMode,
                                       collapse_intron_flank,
                                       gg_theme) {
  time_before <- Sys.time()
  res <- annotation_track_allsamples(df, id, display_region, cds_annotation,
                                     tx_annotation, viewMode, collapse_intron_flank,
                                     gg_theme)

  timer_done_nice_print("-- Mega browser annotation plot done: ", time_before)
  return(res)
}

add_alpha <- function(cols, alpha = 0.7) {
  rgb <- grDevices::col2rgb(cols)
  sprintf("rgba(%d,%d,%d,%.3f)",
          rgb[1, ], rgb[2, ], rgb[3, ], alpha)
}

profile_plotly_gl <- function(dt, frame_colors = frame_color_themes("R"),
                              bar_px = 6, alpha = 0.7) {
  # dt needs: position, count, frame
  stopifnot(is.data.table(dt) || is.data.frame(dt))
  dt <- as.data.table(dt)


  # ensure deterministic order
  dt[, frame := as.character(frame)]
  frame_colors <- add_alpha(frame_colors, alpha)

  p <- plot_ly()
  for (fr in unique(dt$frame)) {
    d <- dt[frame == fr]
    if (nrow(d) == 0L) next

    tt <- paste0(
      "frame: ", fr,
      "<br>pos: ", d$position,
      "<br>count: ", d$count
    )
    # Build vertical segments: (x,0)->(x,count), separated by NA
    xseg <- c(rbind(d$position, d$position, NA_real_))
    yseg <- c(rbind(rep(0, nrow(d)), d$count,   NA_real_))
    textseg <- c(rbind(tt, tt, NA_character_))

    col <- frame_colors[as.integer(fr) + 1L]
    if (is.null(col)) col <- "grey50"

    p <- p %>%
      add_trace(
        x = xseg, y = yseg,
        type = "scattergl", mode = "lines",
        name = fr,
        line = list(color = col, width = bar_px),
        hoverinfo = "text",
        text = textseg,
        showlegend = TRUE
      )
  }

  p %>%
    layout(
      showlegend = FALSE,        # remove legend globally
      margin = margin_megabrowser(),
      paper_bgcolor = "rgba(0,0,0,0)",
      plot_bgcolor  = "rgba(0,0,0,0)",
      xaxis = list(
        showticklabels = FALSE,
        ticks = "",
        showgrid = FALSE,
        zeroline = FALSE,
        title = NULL,
        fixedrange = TRUE
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

summary_track_allsamples <- function(mat, summary_track_type = "area", as_plotly = FALSE) {
  time_before <- Sys.time()
  summary_plot <- profile_plotly_gl(mat)
  timer_done_nice_print("-- Mega browser summary track done: ", time_before)
  return(summary_plot)
}

annotation_track_allsamples <- function(df, id, display_range, annotation,
                                        tx_annotation, viewMode, collapse_intron_flank,
                                        gg_theme) {
  browser()
  gene_model_panel <- createGeneModelPanel(display_range, annotation,
                                           tx_annotation = tx_annotation,
                                           custom_regions = NULL,
                                           viewMode = viewMode, collapse_intron_flank)

  # return(geneModelPanelPlot(gene_model_panel[[1]], attr(gg_theme, "gg_template")) +
  #          theme(plot.margin = margin(l = 23, r = 5)))
  return(geneModelPanelPlotly(gene_model_panel[[1]]))
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
  if (!(enrichment_term %in% valid_enrichment_terms)) {
    stop("Enrichment term is not valid, valid options:\n",
         paste(valid_enrichment_terms, collapse = ", "))
  }
  return(invisible(NULL))
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
