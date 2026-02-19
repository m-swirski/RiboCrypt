allsamples_metadata_clustering <- function(values, plot, enrichment_test_on = "Cluster",
                                           numeric_bins = 5) {

  print("Starting metabrowser clustering info")
  at_least_2_values <- length(unique(values)) > 1
  if (!at_least_2_values) message("Single value analysis not possible, skipping!")

  time_before <- Sys.time()
  req(at_least_2_values)
  row_orders <- row_order(plot)
  orders <- unlist(row_orders, use.names = FALSE)
  clustering_was_done <- length(row_orders) > 1

  if (!clustering_was_done) {
    orders <- order(values) # Order by variable instead of cluster
  }
  # TODO Fix ordering if clustering
  meta <- data.table(grouping = values, order = orders, attr(values, "other_columns"))
  meta <- meta[meta$order, ]
  if (clustering_was_done) {
    meta[, cluster := rep(seq(length(row_orders)), lengths(row_orders))]
  }
  meta[, index := .I]
  attr(meta, "ylab") <- "Enrichment"
  if (enrichment_test_on %in% c("Ratio bins", "Other gene tpm bins")) {
    meta[, cluster := cut(grouping, breaks=numeric_bins)]
    meta[, grouping_numeric_bins_temp := cluster]
    other_cols <- attr(values, "other_columns")
    meta[, grouping := other_cols[,1][[1]]]
    attr(meta, "xlab") <- attr(values, "xlab")
  } else if (is.numeric(meta$grouping)) { #Make numeric bins
      meta[, grouping_numeric_bins := cut(grouping, breaks=numeric_bins)]
    attr(meta, "xlab") <- paste(attr(values, "xlab"), "(Numeric Bins)")
  } else if (enrichment_test_on == "Clusters") {
    attr(meta, "xlab") <- paste(attr(values, "xlab"), "Clusters (K-means)")
  } else {
    attr(meta, "xlab") <- paste0(attr(values, "xlab"), ifelse(clustering_was_done, " Clusters (K-means)", ""))
    attr(meta, "ylab") <- ifelse(clustering_was_done, "Enrichment", "Counts")
  }
  attr(meta, "runIDs") <- attr(values, "runIDs")[orders,]

  enrich_dt <- allsamples_meta_stats(meta)
  timer_done_nice_print("-- metabrowser clustering info done: ", time_before)
  return(list(meta = meta, enrich_dt = enrich_dt))
}

allsamples_sidebar <- function(meta) {
  time_before <- Sys.time()
  gg <- allsamples_sidebar_ggproto(meta)
  columns_to_drop <- c("order", "index", "cluster", attr(meta, "xlab"))
  other_columns <- meta[, !(colnames(meta) %in% columns_to_drop), with = FALSE]
  plot_list <- lapply(other_columns, function(column) {
    if (is.numeric(column)) {
      ggplotly(gg + geom_line(aes(y = column, x = rev(index), fill = NULL)) +
                 coord_flip(), tooltip="text")
    } else {
      ggplotly(gg + geom_raster(aes(fill = column)), tooltip="text")
    }
  })

  res <- subplot(plot_list, nrows = 1) %>% plotly::config(displayModeBar = FALSE)
  timer_done_nice_print("-- metabrowser sidebar done: ", time_before)
  return(res)
}


allsamples_sidebar_plotly <- function(meta) {
  time_before <- Sys.time()

  stopifnot(is.data.table(meta))
  stopifnot(all(c("index", "grouping") %in% names(meta)))
  columns_to_drop <- c("order", "index", "cluster", attr(meta, "xlab"))
  cols <- setdiff(names(meta), columns_to_drop)
  y <- meta[["index"]]

  plot_list <- vector("list", length(cols))
  names(plot_list) <- cols
  for (j in seq_along(cols)) {
    col <- cols[j]
    xj  <- meta[[col]]
    hover_text <- as.character(xj)  # hover text only
    if (is.numeric(xj)) {
      # thin track line; hover shows grouping
      plot_list[[j]] <- plot_ly(
        x = xj, y = y,
        type = "scatter", mode = "lines",
        hoverinfo = "text", text = hover_text,
        line = list(width = 1), showlegend = FALSE
      ) %>%
        layout(
          margin = list(l = 0, r = 0, t = 0, b = 0),
          paper_bgcolor = "rgba(0,0,0,0)",
          plot_bgcolor  = "rgba(0,0,0,0)",
          xaxis = list(showticklabels = FALSE, ticks = "", showgrid = FALSE, zeroline = FALSE, title = NULL,
                       rangeslider = list(visible = FALSE)),
          yaxis = list(showticklabels = FALSE, ticks = "", showgrid = FALSE, zeroline = FALSE, title = NULL)
        )
    } else {
      # categorical track as 1-col heatmap; hover shows grouping
      f <- factor(xj, exclude = NULL)
      z <- matrix(as.integer(f), ncol = 1)

      plot_list[[j]] <- plot_ly(
        x = 1, y = y, z = z,
        type = "heatmap",
        showscale = FALSE,
        hoverinfo = "text",
        text = matrix(hover_text, ncol = 1)
      ) %>%
        layout(
          margin = list(l = 0, r = 0, t = 0, b = 0),
          paper_bgcolor = "rgba(0,0,0,0)",
          plot_bgcolor  = "rgba(0,0,0,0)",
          xaxis = list(showticklabels = FALSE, ticks = "", showgrid = FALSE, zeroline = FALSE, title = NULL,
                       rangeslider = list(visible = FALSE)),
          yaxis = list(showticklabels = FALSE, ticks = "", showgrid = FALSE, zeroline = FALSE, title = NULL)
        )
    }
  }
  res <- subplot(plot_list, nrows = 1, shareY = TRUE, titleX = FALSE, titleY = FALSE) %>%
    plotly::layout(
      xaxis = list(showticklabels = FALSE, ticks = "", showgrid = FALSE, zeroline = FALSE, title = NULL),
      yaxis = list(showticklabels = FALSE, ticks = "", showgrid = FALSE, zeroline = FALSE, title = NULL)
    ) %>%
    plotly::config(displayModeBar = FALSE)

  did_clustering <- length(unique(meta$cluster)) > 1
  if (did_clustering) {
    centers <- meta[, .(
      y_center = (min(index) + max(index)) / 2
    ), by = cluster][order(cluster)]
    centers[, cluster_idx := seq_len(.N)]

    cluster_ann <- lapply(seq_len(nrow(centers)), function(i) {
      list(
        xref = "paper",   # paper coords for x so it stays at the left edge
        yref = "y",       # data coords for y so it matches heatmap rows
        x = -0.02,        # move left of plotting area; tweak as needed
        y = centers$y_center[i],
        text = as.character(centers$cluster_idx[i]),
        hovertext = as.character(centers$cluster[i]),
        captureevents = TRUE,
        showarrow = FALSE,
        textangle = -90,  # sideways
        xanchor = "right",
        yanchor = "middle",
        font = list(size = 12, color = "black")
      )
    })
    res <- res %>% layout(margin = list(l = 30), annotations = cluster_ann)
  }
  timer_done_nice_print("-- metabrowser sidebar done: ", time_before)
  res %>% plotly::layout(margin = list(t = 8, b = 8), yaxis = list(
        autorange = "reversed",
        showticklabels = FALSE,
        ticks = "",
        showgrid = FALSE,
        zeroline = FALSE,
        title = NULL
      )
    )
}


allsamples_sidebar_ggproto <- function(meta) {
  ggplot(meta, aes(y = rev(index), x = factor(1), fill = grouping)) +
    theme_void() +
    labs(x = NULL, y = NULL, title = NULL) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(panel.background = element_rect(fill = "transparent", colour = NA),
          plot.background = element_rect(fill = "transparent", colour = NA),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          plot.margin = unit(c(0, 0, 0, 0), "cm"),
          legend.position = "none")
}


allsamples_meta_stats_shiny <- function(dt) {
  # Add Chi squared significane coloring
  datatable(round(dt, 2)) %>% formatStyle(columns = seq(ncol(dt)),
                                          backgroundColor = styleInterval(c(-3, 3), c('yellow', 'white', 'yellow')))
}

allsamples_meta_stats <- function(meta, attr_xlab = attr(meta, "xlab"), attr_ylab = attr(meta, "ylab")) {
  time_before <- Sys.time()
  print("Starting metabrowser statistics")
  res <- copy(meta)
  res[, index := NULL]
  if ("grouping_numeric_bins" %in% colnames(res)) {
    res[, grouping := grouping_numeric_bins]
  }
  res[, grouping := as.character(grouping)]

  clustering_was_done <- !is.null(res$cluster)
  if (clustering_was_done) {
    concat_table <- table(res$grouping, res$cluster)
    chi_test <- chisq.test(concat_table)
    res <- chi_test$stdres
    tooltipe <- "Chi-squared-stdres: "
  } else {
    concat_table <- table(res$grouping)
    res <- matrix(concat_table, ncol = 1)
    rownames(res) <- names(concat_table)
    colnames(res) <- "counts"
    tooltipe <- "Counts: "
  }
  res <- as.data.frame.matrix(res)
  attr(res, "xlab") <- attr_xlab
  attr(res, "ylab") <- attr_ylab
  attr(res, "tooltip") <- tooltipe

  timer_done_nice_print("-- metabrowser statistics done: ", time_before)
  return(res)
}

allsamples_enrich_bar_plotly <- function(enrich) {
  stopifnot(!is.null(enrich))
  time_before <- Sys.time()

  # Build long table (same semantics as your ggplot version)
  enrich_dt <- data.table::as.data.table(enrich, keep.rownames = TRUE)
  enrich_dt <- suppressWarnings(data.table::melt(enrich_dt, id.vars = "rn"))
  enrich_dt[, variable := factor(as.character(variable))]
  enrich_dt <- enrich_dt[rn != "", ]

  did_enrichment_test <- all(as.character(enrich_dt$variable) != "counts")
  # Tooltips (match your HTML line breaks)
  tip_prefix <- attr(enrich, "tooltip")
  if (is.null(tip_prefix)) tip_prefix <- ""

  xlab_attr <- attr(enrich, "xlab")
  if (is.null(xlab_attr)) xlab_attr <- "Category"

  ylab_attr <- attr(enrich, "ylab")
  if (is.null(ylab_attr)) ylab_attr <- "Value"

  enrich_dt[, tooltip_text := paste0(
    "Category : ", rn, "<br>",
    tip_prefix, round(value, 2), "<br>",
    if (did_enrichment_test) {
      paste0(sub("s \\(K-means\\)", "", xlab_attr), ": ", as.character(variable))
    } else {
      ""
    }
  )]

  # Direct plotly grouped bar chart
  p <- plotly::plot_ly(
    data = enrich_dt,
    x = ~rn,
    y = ~value,
    color = ~variable,
    colors = "Set2",
    type = "bar",
    text = ~tooltip_text,
    hoverinfo = "text"
  ) |>
    plotly::layout(
      barmode = "group",
      xaxis = list(
        title = xlab_attr,
        tickangle = 45,
        titlefont = list(size = 30),
        tickfont  = list(size = 20)
      ),
      yaxis = list(
        title = ylab_attr,
        titlefont = list(size = 30),
        tickfont  = list(size = 20)
      ),
      legend = list(
        title = list(text = "Cluster", font = list(size = 22)),
        font  = list(size = 15)
      )
    )

  # Optional dashed red threshold lines at y = 3 and y = -3
  if (did_enrichment_test) {
    p <- p |>
      plotly::layout(
        shapes = list(
          list(
            type = "line", xref = "paper", yref = "y",
            x0 = 0, x1 = 1, y0 = 3, y1 = 3,
            line = list(color = "red", width = 1, dash = "dash"),
            opacity = 0.7
          ),
          list(
            type = "line", xref = "paper", yref = "y",
            x0 = 0, x1 = 1, y0 = -3, y1 = -3,
            line = list(color = "red", width = 1, dash = "dash"),
            opacity = 0.7
          )
        )
      )
  }
  timer_done_nice_print("-- metabrowser enrichment plot done: ", time_before)
  p
}


allsamples_meta_table <- function(meta_and_clusters) {
  cbind(attr(meta_and_clusters$meta, "runIDs"),
        meta_and_clusters$meta)[, c("index", "order") := NULL]
}
