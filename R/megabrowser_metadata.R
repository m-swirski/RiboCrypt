allsamples_metadata_clustering <- function(values, plot, enrichment_test_on = "Cluster",
                                           numeric_bins = 5) {
  time_before <- Sys.time()
  print("Starting metabrowser clustering info")
  at_least_2_values <- length(unique(values)) > 1
  if (!at_least_2_values) message("Single value analysis not possible, skipping!")

  req(at_least_2_values)
  pdf(NULL) # TODO: Make a better fix for blank pdf write
  row_orders <- suppressWarnings(ComplexHeatmap::row_order(plot))
  dev.off()
  orders <- unlist(row_orders, use.names = FALSE)
  clustering_was_done <- is.list(row_orders)

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
  if (enrichment_test_on %in% c("Ratio bins", "Other gene tpm bins")) {
    meta[, cluster := cut(grouping, breaks=numeric_bins)]
    meta[, grouping_numeric_bins_temp := cluster]
    other_cols <- attr(values, "other_columns")
    meta[, grouping := other_cols[,1][[1]]]
  } else if (is.numeric(meta$grouping)) { #Make numeric bins
      meta[, grouping_numeric_bins := cut(grouping, breaks=numeric_bins)]
    }
  enrich_dt <- allsamples_meta_stats(meta)
  cat("metabrowser clustering info done"); print(round(Sys.time() - time_before, 2))
  return(list(meta = meta, enrich_dt = enrich_dt))
}

allsamples_sidebar <- function(meta) {

  gg <- allsamples_sidebar_ggproto(meta)
  columns_to_drop <- c("order", "index", "cluster")
  other_columns <- meta[, !(colnames(meta) %in% columns_to_drop), with = FALSE]
  plot_list <- lapply(other_columns, function(column) {
    if (is.numeric(column)) {
      ggplotly(gg + geom_line(aes(y = column, x = rev(index), fill = NULL)) +
                 coord_flip(), tooltip="text")
    } else {
      ggplotly(gg + geom_raster(aes(fill = column)), tooltip="text")
    }
  })

  res <- subplot(plot_list, nrows = 1)
  return(res %>% plotly::config(displayModeBar = FALSE) %>%
           layout(margin = list(autoexpand = FALSE, t = 4))
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

allsamples_meta_stats <- function(meta) {
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
  } else {
    concat_table <- table(res$grouping)
    res <- matrix(concat_table, ncol = 1)
    rownames(res) <- names(concat_table)
    colnames(res) <- "counts"
  }
  cat("metabrowser statistics done"); print(round(Sys.time() - time_before, 2))
  return(as.data.frame.matrix(res))
}

allsamples_enrich_bar_plot <- function(enrich) {
  enrich_dt <- as.data.table(enrich, keep.rownames = TRUE)
  enrich_dt <- suppressWarnings(melt(enrich_dt))
  enrich_dt[, variable := factor(as.character(variable))]
  enrich_dt <- enrich_dt[rn != "",]

  enrichment_plot <- ggplot(enrich_dt) +
    geom_bar(aes(x = rn, y = value, fill = variable), stat="identity", position=position_dodge()) +
    theme_minimal() + labs(fill = "Cluster") + xlab("Tissue") + ylab("Enrichment") +
    geom_hline(yintercept = c(3, -3), linetype="dashed", color = "red", linewidth=1) +
    theme(axis.title = element_text(size = 32),
          axis.text.x = element_text(size = (22), angle = 45),
          axis.text.y = element_text(size = (22)),
          legend.text = element_text(size = 22),
          legend.title = element_text(size = 32))
  return(enrichment_plot)
}
