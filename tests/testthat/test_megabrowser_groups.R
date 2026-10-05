megabrowser_groups_fixture <- function() {
  matrix <- matrix(seq_len(24), 4, dimnames = list(NULL, paste0("r", 1:6)))
  attr(matrix, "row_order_list") <- list(c(2L, 5L, 1L), c(6L, 3L, 4L))
  attr(matrix, "km") <- list(cluster = c(1, 1, 2, 2, 1, 2))
  attr(matrix, "ratio") <- 3L
  attr(matrix, "summary_cov") <- data.table::data.table(position = 1:4, count = rowSums(matrix))
  metadata <- data.table::data.table(Run = colnames(matrix), BioProject = "study",
    TISSUE = c("a", "b", "a", "b", "a", "b"), CELL_LINE = c("X", "X", "Y", "Y", "X", "Y"), YEAR = c(2020, 2020, 2021, 2021, NA, 2022))
  values <- stats::setNames(metadata$TISSUE, metadata$Run)
  attr(values, "xlab") <- "TISSUE"
  attr(values, "other_columns") <- metadata[, .(Run, CELL_LINE)]
  attr(values, "runIDs") <- metadata[, .(Run, BioProject)]
  list(table = list(table = matrix, metadata_field = values), metadata = metadata)
}

test_that("sample membership and enrichment follow full-matrix clusters", {
  fixture <- megabrowser_groups_fixture()
  grouped <- allsamples_metadata_clustering(fixture$table, "Clusters")
  expect_equal(grouped$meta$Run, c("r4", "r3", "r6", "r1", "r5", "r2"))
  expect_equal(grouped$meta$cluster, c(2, 2, 2, 1, 1, 1))
  selected <- megabrowser_metadata_enrichment(grouped, fixture$metadata, "CELL_LINE")
  expect_equal(as.integer(attr(selected$enrich_dt, "concat_table")), c(3L, 0L, 0L, 3L))
  expect_equal(selected$meta$Run, grouped$meta$Run)
  expect_equal(selected$meta$cluster, grouped$meta$cluster)
  expect_equal(selected$meta$grouping, c("Y", "Y", "Y", "X", "X", "X"))
  before <- data.table::copy(grouped$meta)
  allsamples_meta_table(grouped)
  expect_identical(grouped$meta, before)
})

test_that("display metadata defers statistics without changing later enrichment", {
  fixture <- megabrowser_groups_fixture()
  original <- allsamples_metadata_clustering(fixture$table, "Clusters")
  calls <- 0L
  stats <- allsamples_meta_stats
  testthat::local_mocked_bindings(allsamples_meta_stats = function(...) {
    calls <<- calls + 1L
    stats(...)
  }, .package = "RiboCrypt")
  deferred <- allsamples_metadata_clustering(fixture$table, "Clusters", compute_stats = FALSE)
  expect_null(deferred$enrich_dt)
  expect_equal(calls, 0L)
  expect_identical(deferred$meta, original$meta)
  expect_equal(megabrowser_metadata_enrichment(deferred, fixture$metadata)$enrich_dt, original$enrich_dt)
  expect_equal(calls, 1L)
  expect_equal(megabrowser_metadata_enrichment(deferred, fixture$metadata, "CELL_LINE"),
               megabrowser_metadata_enrichment(original, fixture$metadata, "CELL_LINE"))
  expect_null(deferred$enrich_dt)
})

test_that("group summaries and memberships use all libraries in the visible original groups", {
  fixture <- megabrowser_groups_fixture()
  grouped <- allsamples_metadata_clustering(fixture$table, "Clusters", compute_stats = FALSE)
  before <- copy(grouped$meta)
  full <- megabrowser_group_membership(fixture$table$table, grouped$meta)
  expect_equal(full$Run, c("r1", "r2", "r5", "r3", "r4", "r6"))
  expect_equal(full$Group, rep(c("1", "2"), each = 3L))
  expect_false(any(c("index", "order", "cluster") %in% names(full)))
  focused <- megabrowser_group_membership(fixture$table$table, grouped$meta, "2")
  expect_equal(focused$Run, c("r3", "r4", "r6"))
  summary <- megabrowser_group_summary(focused)
  expect_true(all(summary$Group == "2" & summary$Libraries == 3L))
  expect_equal(summary[Field == "CELL_LINE", Value], "Y")
  expect_equal(summary[Field == "CELL_LINE", `Agreement (%)`], 100)
  expect_equal(summary[Field == "grouping", `Agreement (%)`], 200 / 3)
  expect_identical(grouped$meta, before)
  expect_equal(megabrowser_metadata_summary(c("X", "Y", NA, ""))$`Agreement (%)`, 25)
  expect_equal(megabrowser_metadata_summary(c("X", "Y", NA, ""))$Missing, 2L)
  expect_equal(megabrowser_group_value(c("", "", NA, "X")), "X")
  expect_equal(megabrowser_group_value(c("", NA)), "(Missing)")
  expect_equal(megabrowser_group_value(c("Y", "X")), "X")
  expect_equal(megabrowser_metadata_summary(c(1, 3, NA))$Value, "2")
  expect_equal(megabrowser_metadata_summary(c(1, 3, NA))$Minimum, 1)
  expect_equal(megabrowser_metadata_summary(c(1, 3, NA))$Maximum, 3)
  expect_true(is.na(megabrowser_metadata_summary(c(NA_real_, NA_real_))$Value))
})

test_that("heatmap preparation preserves exact saved order for matrix and data.table inputs", {
  fixture <- megabrowser_groups_fixture()
  table <- fixture$table$table
  for (input in list(table, set_collection_user_attributes(as.data.table(table), collection_user_attributes(table)))) {
    data <- megabrowser_plotly_heatmap_data(input)
    expected <- megabrowser_ordered_matrix(input)
    expect_identical(data$z, expected[rev(seq_len(nrow(expected))), , drop = FALSE])
    expect_equal(data$x, c(1, 4, 7, 10))
    expect_equal(data$y_range, c(0.5, 6.5))
  }
})

test_that("group CSV exports retain membership and numeric precision", {
  fixture <- megabrowser_groups_fixture()
  grouped <- allsamples_metadata_clustering(fixture$table, "Clusters", compute_stats = FALSE)
  membership <- megabrowser_group_membership(fixture$table$table, grouped$meta, "2")
  membership[, score := c(1 / 3, 2 / 3, NA_real_)]
  summary <- megabrowser_group_summary(membership)
  testthat::local_mocked_bindings(downloadHandler = function(filename, content, contentType) {
    list(filename = filename, content = content, type = contentType)
  }, .package = "RiboCrypt")
  for (name in c("group-summary", "library-membership")) {
    expected <- copy(if (name == "group-summary") summary else membership)
    data.table::setindexv(expected, NULL)
    handler <- megabrowser_csv_download(function() expected, name)
    file <- tempfile(fileext = ".csv")
    handler$content(file)
    expect_equal(handler$filename(), paste0("megabrowser-", name, ".csv"))
    expect_equal(handler$type, "text/csv")
    downloaded <- fread(file, colClasses = c(Group = "character"))
    expect_equal(lapply(downloaded, as.vector), lapply(expected, as.vector))
    unlink(file)
  }
})

test_that("collapsed matrices contain means of existing groups and preserve the full state", {
  fixture <- megabrowser_groups_fixture()
  original <- unserialize(serialize(fixture$table, NULL))
  grouped <- allsamples_metadata_clustering(fixture$table, "Clusters")
  display <- megabrowser_collapsed_display(fixture$table$table, grouped$meta)
  expect_equal(unname(display$table[, 1]), rowMeans(fixture$table$table[, c(2, 5, 1)]))
  expect_equal(unname(display$table[, 2]), rowMeans(fixture$table$table[, c(6, 3, 4)]))
  expect_equal(display$meta$cluster, c("2", "1"))
  expect_equal(display$meta$libraries, c(3L, 3L))
  expect_equal(attr(display$table, "summary_cov"), attr(fixture$table$table, "summary_cov"))
  expect_equal(attr(display$table, "ratio"), 3L)
  expect_equal(fixture$table, original)
  heatmap <- megabrowser_plotly_heatmap_data(display$table)
  expect_equal(heatmap$y, 1:2)
  expect_equal(heatmap$y_range, c(0.5, 2.5))
  expect_equal(heatmap$type, "heatmap")
  plot <- plotly::plotly_build(megabrowser_plotly_heatmap(display$table, c("white", "blue"), megabrowserHeatmapPlotlyTemplate()))
  expect_equal(plot$x$data[[1]]$type, "heatmap")
  expect_equal(plot$x$data[[1]]$y0, 1)
  expect_equal(plot$x$layout$yaxis$dtick, 1)
  expect_equal(as.numeric(heatmap$z[1, ]), rowMeans(fixture$table$table[, c(2, 5, 1)]))
  static <- megabrowser_complex_heatmap(display$table, c("white", "blue"))
  expect_equal(dim(static@matrix), c(2L, 4L))
  expect_equal(as.numeric(static@matrix[1, ]), as.numeric(heatmap$z[nrow(heatmap$z), ]))
})

test_that("numeric bins, categories, missing metadata and single-row matrices collapse safely", {
  fixture <- megabrowser_groups_fixture()
  bins <- megabrowser_numeric_bins(c(1, 1, NA, Inf))
  expect_equal(as.character(bins), c("All", "All", "(Missing)", "(Missing)"))
  attr(fixture$table$table, "row_order_list") <- list(1:6)
  grouped <- allsamples_metadata_clustering(fixture$table, "TISSUE")
  display <- megabrowser_collapsed_display(fixture$table$table, grouped$meta)
  expect_equal(ncol(display$table), 2L)
  expect_equal(sum(display$meta$libraries), 6L)
  numeric <- megabrowser_metadata_enrichment(grouped, fixture$metadata, "YEAR")
  expect_equal(sum(attr(numeric$enrich_dt, "concat_table")), 6L)
  old_attributes <- attributes(fixture$table$metadata_field)
  fixture$table$metadata_field <- seq_len(6)
  attributes(fixture$table$metadata_field) <- old_attributes
  attr(fixture$table$metadata_field, "other_columns") <- fixture$metadata[, .(Run, TISSUE)]
  grouped <- allsamples_metadata_clustering(fixture$table, "Ratio bins")
  expect_equal(sum(attr(grouped$enrich_dt, "concat_table")), 6L)
  display <- megabrowser_collapsed_display(fixture$table$table[1, , drop = FALSE], grouped$meta)
  expect_equal(nrow(display$table), 1L)
  expect_equal(ncol(display$table), length(unique(grouped$meta$cluster)))
})

test_that("focus remaps display rows without changing clusters or sample statistics", {
  fixture <- megabrowser_groups_fixture()
  grouped <- allsamples_metadata_clustering(fixture$table, "Clusters")
  original <- unserialize(serialize(fixture$table, NULL))
  focused <- megabrowser_focused_display(fixture$table$table, grouped$meta, "2")
  expect_equal(colnames(focused$table), c("r3", "r4", "r6"))
  expect_equal(attr(focused$table, "row_order_list"), list(c(3L, 1L, 2L)))
  expect_equal(attr(focused$table, "km")$cluster, c(2, 2, 2))
  expect_equal(focused$meta$Run, c("r4", "r3", "r6"))
  expect_equal(focused$meta$index, 1:3)
  expect_equal(attr(focused$table, "summary_cov"), attr(fixture$table$table, "summary_cov"))
  expect_equal(fixture$table, original)
  collapsed <- megabrowser_collapsed_display(focused$table, focused$meta)
  expect_equal(as.numeric(collapsed$table), rowMeans(fixture$table$table[, c(6, 3, 4)]))
  expect_equal(collapsed$meta$libraries, 3L)
  sidebar <- plotly::plotly_build(allsamples_sidebar_plotly(collapsed$meta))
  expect_equal(sidebar$x$layout$annotations[[1]]$text, "2")
  expect_equal(megabrowser_display_counts(fixture$table$table, focused), "3 of 6 libraries | 3 rows | 4 positions")
  expect_equal(megabrowser_display_counts(fixture$table$table, collapsed), "3 of 6 libraries | 1 row | 4 positions")
  expect_identical(megabrowser_focused_display(fixture$table$table, grouped$meta, "missing")$table, fixture$table$table)
  expect_identical(megabrowser_focused_display(fixture$table$table, grouped$meta, c("2", "1"))$table, fixture$table$table)
})

test_that("moderate matrices use canvas while large matrices retain WebGL", {
  expect_identical(megabrowser_heatmap_renderer(2e6), "heatmap")
  expect_identical(megabrowser_heatmap_renderer(2e6 + 1), "heatmapgl")
  expect_identical(megabrowser_heatmap_renderer(2e6 + 1, collapsed = TRUE), "heatmap")
  matrix <- matrix(0, 1001, 2000)
  attr(matrix, "km") <- list(cluster = rep(1, 2000))
  attr(matrix, "row_order_list") <- list(1:2000)
  expect_identical(megabrowser_plotly_heatmap_data(matrix)$type, "heatmapgl")
  expect_false("dtick" %in% names(megabrowser_heatmap_layout(megabrowser_plotly_heatmap_data(matrix))$yaxis))
})

test_that("full and collapsed heatmaps hide x-axis decorations without changing ranges", {
  fixture <- megabrowser_groups_fixture()
  grouped <- allsamples_metadata_clustering(fixture$table, "Clusters")
  collapsed <- megabrowser_collapsed_display(fixture$table$table, grouped$meta)
  for (table in list(fixture$table$table, collapsed$table)) {
    for (template in list(NULL, megabrowserHeatmapPlotlyTemplate())) {
      plot <- plotly::plotly_build(megabrowser_plotly_heatmap(table, c("white", "red"), template))
      axis <- plot$x$layout$xaxis
      expect_false(axis$showticklabels)
      expect_identical(axis$ticks, "")
      expect_false(axis$showgrid)
      expect_false(axis$zeroline)
      expect_false(axis$showline)
      expect_false(axis$automargin)
      expect_equal(axis$range, megabrowser_full_x_range(table = table))
      expect_equal(plot$x$layout$margin$b, 8)
    }
  }
})

test_that("sidebar reset and fractional zoom preserve exact heatmap boundaries", {
  for (rows in c(1L, 2L, 5L, 1754L)) {
    expected <- c(rows + 0.5, 0.5)
    reset <- mb_y_relayout_from_event(list(r = c(0.5, rows + 0.5)), y_max = rows)
    expect_equal(reset$yaxis$range, expected)
    expect_false(reset$yaxis$autorange)
    expect_equal(mb_y_relayout_from_event(list(auto = TRUE), y_max = rows)$yaxis$range, expected)
  }
  expect_equal(mb_sidebar_y_range(1.25, 3.75, y_max = 5), c(4.75, 2.25))
})

test_that("sidebar row boundaries match full, collapsed and single-group heatmaps", {
  fixture <- megabrowser_groups_fixture()
  grouped <- allsamples_metadata_clustering(fixture$table, "Clusters")
  focused <- megabrowser_focused_display(fixture$table$table, grouped$meta, "2")
  displays <- list(list(table = fixture$table$table, meta = grouped$meta),
    megabrowser_collapsed_display(fixture$table$table, grouped$meta),
    megabrowser_collapsed_display(focused$table, focused$meta))
  for (display in displays) {
    sidebar <- plotly::plotly_build(allsamples_sidebar_plotly(display$meta))
    heatmap <- megabrowser_plotly_heatmap_data(display$table)
    expect_false(sidebar$x$layout$yaxis$autorange)
    expect_equal(sidebar$x$layout$yaxis$range, rev(heatmap$y_range))
    expect_equal(sidebar$x$layout$margin$t, megabrowser_heatmap_margins()$t)
    expect_equal(sidebar$x$layout$margin$b, megabrowser_heatmap_margins()$b)
  }
})

test_that("nonadjacent group focus preserves saved order and original group labels", {
  fixture <- megabrowser_groups_fixture()
  attr(fixture$table$table, "row_order_list") <- list(c(2L, 5L), c(1L, 6L), c(3L, 4L))
  attr(fixture$table$table, "km")$cluster <- c(2L, 1L, 3L, 3L, 1L, 2L)
  grouped <- allsamples_metadata_clustering(fixture$table, "Clusters")
  focused <- megabrowser_focused_display(fixture$table$table, grouped$meta, c("3", "1"))
  expect_equal(colnames(focused$table), c("r2", "r3", "r4", "r5"))
  expect_equal(attr(focused$table, "row_order_list"), list(c(1L, 4L), c(2L, 3L)))
  expect_equal(attr(focused$table, "km")$cluster, c(1L, 3L, 3L, 1L))
  collapsed <- megabrowser_collapsed_display(focused$table, focused$meta)
  expect_equal(colnames(collapsed$table), c("1", "3"))
  expect_equal(collapsed$meta$cluster, c("3", "1"))
  expect_equal(unname(collapsed$table[, 1]), rowMeans(fixture$table$table[, c(2, 5)]))
  expect_equal(unname(collapsed$table[, 2]), rowMeans(fixture$table$table[, c(3, 4)]))
})

test_that("changing enrichment or collapse never loads coverage or runs clustering again", {
  fixture <- megabrowser_groups_fixture()
  computations <- plots <- output_builds <- 0L
  reset_plot <- megabrowser_mid_reset_plot
  controls <- list(table_hash = "groups-test", table_plot_hash = "groups-plot-test",
    enrichment_term = "Clusters", plotType = "plotly", display_region = NULL)
  testthat::local_mocked_bindings(
    allsamples_observer_controller = function(...) NULL,
    mb_controller_shiny = function(...) controls,
    compute_collection_table_shiny = function(...) {computations <<- computations + 1L; fixture$table},
    get_megabrowser_annotation_plot_shiny = function(...) plotly::plot_ly(x = 1, y = 1, type = "scatter", mode = "markers"),
    mb_plot_object_shiny = function(table_obj, ...) {plots <<- plots + 1L; megabrowser_plotly_heatmap(table_obj, c("white", "blue"))},
    megabrowser_mid_reset_plot = function(...) {output_builds <<- output_builds + 1L; reset_plot(...)},
    .package = "RiboCrypt")
  shiny::testServer(browser_allsamp_server, args = list(id = "mb", all_exp = NULL, df = NULL,
    experiments = NULL, gene_name_list = NULL, tx = NULL, cds = NULL, org = NULL,
    motif_name_list = NULL, metadata = fixture$metadata, browser_options = NULL, rv = NULL), {
    session$setInputs(go = 0, plotType = "plotly", collapsed_clusters = FALSE, enrichment_metadata = "grouping")
    session$setInputs(go = 1)
    plot_object()
    original_output <- output$myPlotlyPlot
    expect_true(nzchar(original_output))
    expect_equal(output_builds, 1L)
    initial_stats <- meta_and_clusters()$enrich_dt
    expect_equal(computations, 1L)
    expect_equal(plots, 1L)
    session$setInputs(enrichment_metadata = "CELL_LINE")
    expect_false(identical(initial_stats, meta_and_clusters()$enrich_dt))
    expect_equal(computations, 1L)
    expect_equal(plots, 1L)
    full_stats <- meta_and_clusters()$enrich_dt
    session$setInputs(collapsed_clusters = TRUE)
    expect_equal(ncol(display_table()$table), 2L)
    expect_equal(output$display_counts, "6 libraries | 2 rows | 4 positions")
    plot_object()
    expect_equal(plots, 2L)
    expect_equal(computations, 1L)
    expect_identical(meta_and_clusters()$enrich_dt, full_stats)
    expect_equal(ncol(table()$table), 6L)
    expect_equal(nrow(allsamples_meta_table(meta_and_clusters())), 6L)
    expect_false(identical(output$myPlotlyPlot, original_output))
    expect_equal(output_builds, 2L)
    session$setInputs(collapsed_clusters = FALSE)
    expect_identical(output$myPlotlyPlot, original_output)
    expect_equal(output_builds, 2L)
    expect_equal(plots, 2L)
    expect_equal(computations, 1L)
    session$setInputs(visible_groups = "2")
    expect_equal(ncol(display_table()$table), 3L)
    expect_identical(meta_and_clusters()$enrich_dt, full_stats)
    expect_equal(nrow(allsamples_meta_table(meta_and_clusters())), 6L)
    expect_equal(computations, 1L)
    plot_object()
    expect_equal(plots, 3L)
    session$setInputs(collapsed_clusters = TRUE)
    expect_equal(ncol(display_table()$table), 1L)
    expect_equal(display_table()$meta$libraries, 3L)
    expect_identical(meta_and_clusters()$enrich_dt, full_stats)
    expect_equal(computations, 1L)
    session$setInputs(visible_groups = character(), collapsed_clusters = FALSE)
    expect_identical(output$myPlotlyPlot, original_output)
    expect_equal(computations, 1L)
  })
})

test_that("view controls are scoped and layout preferences are client-only", {
  controls <- as.character(megabrowser_view_controls(shiny::NS("mb")))
  expect_match(controls, 'id="mb-collapsed_clusters"', fixed = TRUE)
  expect_match(controls, 'id="mb-enrichment_metadata"', fixed = TRUE)
  expect_match(controls, 'aria-label="Reset zoom"', fixed = TRUE)
  expect_match(controls, 'class="mega-height-control"', fixed = TRUE)
  layout <- as.character(renderMegabrowser("plotly", shiny::NS("mb"), height = "var(--mega-height, 700px)"))
  expect_match(layout, "calc(var(--mega-height, 700px) * 0.750000)", fixed = TRUE)
  static_layout <- as.character(renderMegabrowser("ggplot2", shiny::NS("mb"), height = "var(--mega-height, 700px)"))
  expect_match(static_layout, "calc(var(--mega-height, 700px) * 1)", fixed = TRUE)
  expect_false(grepl("Shiny.setInputValue", fetchJS("megabrowser_layout.js"), fixed = TRUE))
})
