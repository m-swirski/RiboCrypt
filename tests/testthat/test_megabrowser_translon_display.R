translon_display_fixture <- function() {
  raw <- matrix(100, 10, 3, dimnames = list(NULL, c("a", "b", "c")))
  raw[1:2, ] <- rep(c(4, 8, 12), each = 2)
  raw[c(5:6, 9:10), ] <- rep(c(2, 6, 10), each = 4)
  table <- matrix(99, 5, 3, dimnames = list(NULL, colnames(raw)))
  attr(table, "ratio") <- 2
  attr(table, "km") <- list(cluster = c(1, 1, 2))
  attr(table, "row_order_list") <- list(c(1, 2), 3)
  attr(table, "clusters") <- 2
  attr(table, "summary_cov") <- data.table(position = 1:10)
  regions <- list(R1 = list(labels = "T1, TC1", kind = "Translon", ranges = IRanges::IRanges(1, 2)),
    R2 = list(labels = "CDS", kind = "CDS", ranges = IRanges::IRanges(c(5, 9), c(6, 10))),
    R3 = list(labels = "clean_cds", parent_labels = "CDS", kind = "clean_cds", ranges = IRanges::IRanges(9, 10)),
    R4 = megabrowser_user_region("3", "UR1", 10))
  list(raw = raw, table = table, regions = regions,
    meta = data.table(Run = c("a", "b", "c"), cluster = c(1, 1, 2), grouping = "x"))
}

test_that("translon columns are exact length-normalized raw densities with saved clustering", {
  f <- translon_display_fixture()
  display <- megabrowser_translon_display(f$table, f$raw[, 3:1], f$regions)
  expect_equal(unname(as.vector(display)), c(4, 100, 2, 8, 100, 6, 12, 100, 10))
  expect_equal(dimnames(display), list(c("R1", "R4", "R3"), c("a", "b", "c")))
  expect_identical(attr(display, "row_order_list"), attr(f$table, "row_order_list"))
  expect_equal(attr(display, "translon_summary"), c(R1 = 8, R4 = 100, R3 = 6))
  expect_equal(attr(display, "ratio"), 1)
  expect_equal(as.vector(f$table), rep(99, 15))
  expect_equal(attr(f$table, "ratio"), 2)
  expect_equal(attr(display, "translon_regions")$Length_nt, c(2L, 1L, 2L))
})

test_that("translon and cluster collapse compose without changing the full analysis matrix", {
  f <- translon_display_fixture()
  display <- megabrowser_translon_display(f$table, f$raw, f$regions)
  collapsed <- megabrowser_collapsed_display(display, f$meta)
  groups <- megabrowser_display_groups(f$meta, colnames(display))
  for (group in names(groups)) expect_equal(as.vector(collapsed$table[, group]), unname(rowMeans(display[, groups[[group]], drop = FALSE])))
  expect_equal(ncol(collapsed$table), 2)
  expect_true(attr(collapsed$table, "collapsed_translons"))
  expect_true(attr(collapsed$table, "collapsed_clusters"))
  focused <- megabrowser_focused_display(display, f$meta, "1")
  expect_equal(colnames(focused$table), c("a", "b"))
  expect_equal(nrow(focused$table), 3)
  expect_match(megabrowser_display_counts(f$table, collapsed), "3 libraries.*2 rows.*3 regions")
  scored <- megabrowser_translon_cluster_score(collapsed)
  expect_equal(unname(scored$table["R1", "1"]), log2(6 / 8))
  expect_equal(unname(scored$table["R1", "2"]), log2(12 / 8))
  focused_score <- megabrowser_translon_cluster_score(megabrowser_collapsed_display(focused$table, focused$meta))
  expect_equal(as.vector(focused_score$table[, "1"]), as.vector(scored$table[, "1"]))
  expect_identical(attr(focused_score$table, "translon_summary"), attr(scored$table, "translon_summary"))
})

test_that("equal-width region axes, hover labels and reset geometry work with plot templates", {
  f <- translon_display_fixture()
  display <- megabrowser_translon_display(f$table, f$raw, f$regions)
  data <- megabrowser_plotly_heatmap_data(display)
  expect_equal(data$x, 1:3)
  expect_equal(data$x_range, c(0.5, 3.5))
  expect_equal(megabrowser_full_x_range(GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 500)), display), c(0.5, 3.5))
  expect_equal(data$text[1, ], c("T1, TC1", "UR1", "clean_cds"))
  expect_equal(data$z[1, ], c(R1 = 4, R4 = 100, R3 = 2))
  p <- plotly::plotly_build(megabrowser_plotly_heatmap(display, c("white", "blue"), megabrowserHeatmapPlotlyTemplate()))
  expect_equal(p$x$layout$xaxis$range, c(0.5, 3.5))
  expect_false(p$x$layout$xaxis$showticklabels)
  expect_match(p$x$data[[1]]$hovertemplate, "Mean coverage")
  expect_equal(dim(p$x$data[[1]]$text), c(3, 3))
  top <- plotly::plotly_build(megabrowser_translon_summary_plot(display))
  bottom <- plotly::plotly_build(megabrowser_translon_annotation_plot(display))
  expect_equal(top$x$layout$xaxis$range, c(0.5, 3.5))
  expect_equal(bottom$x$layout$xaxis$range, c(0.5, 3.5))
  expect_equal(as.numeric(top$x$data[[1]]$y), c(8, 100, 6))
  expect_equal(as.character(bottom$x$data[[1]]$text), c("uORF1", "UR1", "clean CDS"))
  expect_null(top$width)
  expect_null(bottom$width)
  expect_equal(bottom$x$layout$margin$l, p$x$layout$margin$l)
  expect_equal(bottom$x$layout$margin$r, p$x$layout$margin$r)
  expect_equal(megabrowser_reset_range_shiny(function() list(display_region = NULL), function() list(table = display)), c(0.5, 3.5))
  expect_s4_class(megabrowser_complex_heatmap(display, c("white", "blue")), "Heatmap")
})

test_that("invalid raw coordinates and coverage give actionable validation", {
  f <- translon_display_fixture()
  expect_error(megabrowser_translon_display(f$table, f$raw[1:5, ], f$regions), "coordinates")
  f$raw[1, 1] <- NA_real_
  expect_error(megabrowser_translon_display(f$table, f$raw, f$regions), "finite and non-negative")
})

test_that("clean CDS replaces only its parent and regions follow transcript order", {
  f <- translon_display_fixture()
  ordered <- megabrowser_ordered_translon_regions(f$regions[c(3, 2, 4, 1)])
  expect_equal(names(ordered), c("R1", "R4", "R3"))
  expect_false(any(vapply(ordered, function(r) r$kind == "CDS", logical(1))))
  expect_equal(names(megabrowser_ordered_translon_regions(f$regions[c(1, 2, 4)])), c("R1", "R4", "R2"))
  f$regions$R3$ranges <- IRanges::IRanges()
  expect_equal(tail(names(megabrowser_ordered_translon_regions(f$regions)), 1), "R3")
  f$regions$R3$parent_start <- 5
  f$regions$R5 <- megabrowser_user_region("10", "UR2", 10)
  expect_equal(names(megabrowser_ordered_translon_regions(f$regions)), c("R1", "R4", "R3", "R5"))
})

test_that("region-relative cluster colours reveal equal fold changes in weak and strong regions", {
  raw <- matrix(c(100, 5, 200, 10, 300, 15), 2, 3, dimnames = list(c("R1", "R2"), c("a", "b", "c")))
  attr(raw, "translon_summary") <- rowMeans(raw)
  attr(raw, "translon_regions") <- data.table(Labels = c("strong", "weak"))
  attr(raw, "collapsed_translons") <- TRUE
  attr(raw, "row_order_list") <- list(1, 2, 3)
  attr(raw, "km") <- list(cluster = 1:3)
  scored <- megabrowser_translon_cluster_score(list(table = raw))$table
  expect_equal(as.vector(scored[1, ]), as.vector(scored[2, ]))
  expect_equal(as.vector(scored[2, ]), log2(c(0.5, 1, 1.5)))
  expect_equal(as.vector(attr(scored, "translon_density")), as.vector(raw))
  expect_equal(attr(scored, "translon_fold_change")[2, 3], 1.5)
  p <- plotly::plotly_build(get_meta_browser_plot(scored, "default (White-Blue)", template = megabrowserHeatmapPlotlyTemplate()))
  expect_equal(p$x$data[[1]]$zmin, -1)
  expect_equal(p$x$data[[1]]$zmax, 1)
  expect_true(p$x$data[[1]]$showscale)
  expect_equal(p$x$data[[1]]$customdata[3, 2], 15)
  expect_match(p$x$data[[1]]$text[3, 2], "1.5")
  attr(raw, "translon_summary") <- c(R1 = 0, R2 = 10)
  raw[2, ] <- c(0, 10, 100)
  zero <- megabrowser_translon_cluster_score(list(table = raw))$table
  expect_true(all(is.na(zero[1, ])))
  expect_equal(as.vector(zero[2, ]), c(-1, 0, 1))
  expect_equal(attr(zero, "translon_fold_change")[2, 1], 0)
  expect_equal(attr(zero, "translon_density")[2, 3], 100)
})

test_that("shared region workspace is lazy, reuses raw coverage and resets custom coordinates", {
  reads <- 0L
  f <- translon_display_fixture()
  panel <- data.table(gene_names = c("T", "CDS"), type = "cds", rect_starts = c(1L, 5L), rect_ends = c(2L, 10L))
  testthat::local_mocked_bindings(load_collection = function(path, columns) {
    force(path)
    force(columns)
    reads <<- reads + 1L
    f$raw
  },
    megabrowser_translon_panel = function(...) panel, .package = "RiboCrypt")
  shiny::testServer(function(input, output, session) {
    controller <- reactive(list(table_hash = input$hash, table_plot_hash = "plot", table_path = "mock", annotation = list(CDS = NULL)))
    workspace <- megabrowser_translon_workspace(controller, reactive(list(table = f$table)))
  }, {
    session$setInputs(hash = "first")
    expect_equal(reads, 0L)
    workspace$raw()
    workspace$regions()
    expect_equal(reads, 1L)
    workspace$custom(list(megabrowser_user_region("3:4", "UR1", 10)))
    expect_equal(tail(workspace$regions(), 1)[[1]]$labels, "UR1")
    workspace$raw()
    expect_equal(reads, 1L)
    session$setInputs(hash = "next")
    expect_length(workspace$custom(), 0)
    expect_equal(reads, 1L)
    workspace$raw()
    expect_equal(reads, 2L)
  })
})

test_that("relative scores respect both selected palettes, colour zoom and legend swatches", {
  f <- translon_display_fixture()
  table <- megabrowser_translon_display(f$table, f$raw, f$regions)
  scored <- megabrowser_translon_cluster_score(megabrowser_collapsed_display(table, f$meta))$table
  plain <- scored
  attr(plain, "translon_score") <- FALSE
  scales <- list()
  for (theme in c("default (White-Blue)", "Matrix (black,green,red)")) for (zoom in c(1, 8)) {
    colors <- megabrowser_heatmap_colors(theme, zoom)
    p <- plotly::plotly_build(get_meta_browser_plot(scored, theme, zoom, template = megabrowserHeatmapPlotlyTemplate()))
    regular <- plotly::plotly_build(get_meta_browser_plot(plain, theme, zoom, template = megabrowserHeatmapPlotlyTemplate()))
    expect_identical(p$x$data[[1]]$colorscale, regular$x$data[[1]]$colorscale)
    expect_equal(p$x$data[[1]]$z, regular$x$data[[1]]$z)
    expect_equal(c(p$x$data[[1]]$zmin, p$x$data[[1]]$zmax), c(-1, 1))
    mapping <- circlize::colorRamp2(seq(-1, 1, length.out = length(colors)), colors, space = "RGB")
    static <- megabrowser_complex_heatmap(scored, colors)
    expect_equal(static@matrix_color_mapping@col_fun(c(-1, 0, 1)), mapping(c(-1, 0, 1)))
    legend <- as.character(megabrowser_score_legend(colors))
    for (color in mapping(c(-1, 0, 1))) expect_match(legend, color, fixed = TRUE)
    scales[[length(scales) + 1L]] <- p$x$data[[1]]$colorscale
  }
  expect_false(identical(scales[[1]], scales[[2]]))
  expect_false(identical(scales[[1]], scales[[3]]))
})
