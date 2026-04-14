test_that("get_ratio_interval parses and validates inputs", {
  expect_equal(RiboCrypt:::get_ratio_interval("5:10;20:30"), c(5, 10, 20, 30))
  expect_equal(RiboCrypt:::get_ratio_interval("7"), c(7, 7))
  expect_error(RiboCrypt:::get_ratio_interval("5:10:15"))
  expect_error(suppressWarnings(RiboCrypt:::get_ratio_interval("a:b")))
})

test_that("validate_enrichment_term uses Shiny validation for invalid inputs", {
  expect_no_error(
    RiboCrypt:::validate_enrichment_term(
      enrichment_term = "TISSUE",
      clusters = 1,
      ratio_interval = NULL,
      other_gene = NULL,
      metadata_field = c("TISSUE", "CELL_LINE")
    )
  )

  err <- tryCatch(
    {
      RiboCrypt:::validate_enrichment_term(
        enrichment_term = "not-valid",
        clusters = 1,
        ratio_interval = NULL,
        other_gene = NULL,
        metadata_field = c("TISSUE", "CELL_LINE")
      )
      NULL
    },
    error = function(e) e
  )

  expect_s3_class(err, "shiny.silent.error")
  expect_match(conditionMessage(err), "Enrichment term is not valid")
  expect_match(conditionMessage(err), "TISSUE")
})

test_that("module controller helpers parse events and observatory kickoff state", {
  session <- list(
    rootScope = function() {
      list(input = list(
        `plotly_click-mb_enrich` = "{\"x\":\"term\",\"customdata\":\"2\"}",
        `plotly_relayout-mb_mid` = list(`xaxis.range[0]` = 10)
      ))
    }
  )

  click_event <- RiboCrypt:::get_plotly_session_event(session, "plotly_click", "mb_enrich")
  expect_equal(click_event$x, "term")
  expect_equal(click_event$customdata, "2")

  relayout_event <- RiboCrypt:::get_plotly_session_event(session, "plotly_relayout", "mb_mid")
  expect_equal(relayout_event[["xaxis.range[0]"]], 10)

  input <- list(gene = "GENE1", tx = "TX1")
  state <- list(view = "browser", browser = list(go = TRUE))
  expect_true(RiboCrypt:::observatory_browser_ready_to_kickoff(state, input, list(A = "SRR1")))
  expect_false(RiboCrypt:::observatory_browser_ready_to_kickoff(state, list(gene = "", tx = "TX1"), list(A = "SRR1")))
  expect_false(RiboCrypt:::observatory_browser_ready_to_kickoff(list(view = "selector", browser = list(go = TRUE)), input, list(A = "SRR1")))
})

test_that("multiSampleBinRows bins rows as expected", {
  mat <- matrix(1:20, nrow = 10, ncol = 2)
  binned <- RiboCrypt:::multiSampleBinRows(mat, ratio = 3)
  expect_equal(nrow(binned), 4)
  expect_equal(binned[1, ], colSums(mat[1:3, ]))

  mat_small <- matrix(1:6, nrow = 3, ncol = 2)
  expect_equal(RiboCrypt:::multiSampleBinRows(mat_small, ratio = 3), mat_small)
})

test_that("summary_track returns expected columns and values", {
  mat <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 3, ncol = 2)
  res <- RiboCrypt:::summary_track(mat)
  expect_true(all(c("count", "position", "frame") %in% names(res)))
  expect_equal(res$count, c(1 + 4, 2 + 5, 3 + 6))
  expect_equal(res$position, 1:3)
})

test_that("megabrowser_full_x_range prefers display range span", {
  table <- data.table::data.table(a = 1:150, b = 151:300)
  summary_cov <- data.table::data.table(
    count = seq_len(1350),
    position = seq_len(1350),
    frame = factor((seq_len(1350) - 1) %% 3)
  )
  data.table::setattr(table, "summary_cov", summary_cov)
  display_range <- GenomicRanges::GRangesList(
    tx = GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 1350), "+")
  )

  expect_equal(RiboCrypt:::megabrowser_full_x_range(display_range, table), c(1, 1350))
})

test_that("margin_megabrowser returns expected margins", {
  m <- RiboCrypt:::margin_megabrowser()
  expect_equal(m$l, 30)
  expect_equal(m$r, 100)
  expect_equal(m$t, 0)
  expect_equal(m$b, 0)
})

test_that("allsamples_meta_stats computes counts without clusters", {
  meta <- data.table::data.table(grouping = c("a", "a", "b"), order = 1:3, index = seq(3))
  attr(meta, "xlab") <- "X"
  attr(meta, "ylab") <- "Y"
  res <- RiboCrypt:::allsamples_meta_stats(meta)
  expect_equal(res["a", "counts"], 2)
  expect_equal(res["b", "counts"], 1)
  expect_equal(attr(res, "xlab"), "X")
  expect_equal(attr(res, "ylab"), "Y")
})

test_that("add_alpha converts colors to rgba with alpha", {
  res <- RiboCrypt:::add_alpha(c("red", "blue"), alpha = 0.5)
  expect_true(all(grepl("^rgba\\(", res)))
  expect_true(all(grepl(",0.500\\)$", res)))
})

test_that("profile_plotly_gl returns a plotly object", {
  dt <- data.table::data.table(
    position = 1:6,
    count = c(1, 2, 3, 4, 5, 6),
    frame = factor(rep(0:2, length.out = 6))
  )
  p <- RiboCrypt:::profile_plotly_gl(dt, bar_px = 2, alpha = 0.6)
  expect_true(inherits(p, "plotly"))
  built <- plotly::plotly_build(p)
  expect_equal(unname(built$x$layout$xaxis$range), c(1, 6))
  expect_identical(built$x$layout$xaxis$autorange, FALSE)
})

test_that("summary_track_allsamples uses browser columns template with zoom switch", {
  dt <- data.table::data.table(
    position = 1:6,
    count = c(1, 2, 3, 4, 5, 6),
    frame = factor(rep(0:2, length.out = 6))
  )
  p <- RiboCrypt:::summary_track_allsamples(
    dt,
    template = RiboCrypt:::covPanelColumnsPlotlyTemplate()
  )

  expect_true(inherits(p, "plotly"))
  expect_identical(p$x$source, "mb_top")
  expect_identical(p$x$config$doubleClick, FALSE)

  built <- plotly::plotly_build(p)
  expect_equal(unname(built$x$layout$xaxis$range), c(1, 6))
  expect_identical(built$x$layout$xaxis$autorange, FALSE)
  expect_equal(unname(built$x$layout$margin$l), 30)
  expect_equal(unname(built$x$layout$margin$r), 100)
  expect_identical(built$x$layout$paper_bgcolor, "rgba(0,0,0,0)")
  expect_identical(built$x$layout$plot_bgcolor, "rgba(0,0,0,0)")
  expect_identical(built$x$layout$hovermode, "closest")
  expect_null(built$x$layout$yaxis$title)
  expect_identical(built$x$layout$yaxis$fixedrange, TRUE)
  expect_identical(built$x$layout$yaxis$showticklabels, FALSE)
  expect_gte(length(built$x$data), 6)
  expect_identical(built$x$data[[1]]$type, "scattergl")
  expect_identical(built$x$data[[1]]$visible, "legendonly")
  expect_identical(built$x$data[[1]]$meta$rc_columns_switch, "gl_subset")
  expect_identical(built$x$data[[4]]$type, "scatter")
  expect_identical(built$x$data[[4]]$meta$rc_columns_switch, "line")
})

test_that("summary_track_allsamples handles empty summary data without warnings", {
  empty_dt <- data.table::data.table(
    count = numeric(),
    position = numeric(),
    frame = factor(levels = c("0", "1", "2"))
  )

  expect_no_warning(
    p <- RiboCrypt:::summary_track_allsamples(
      empty_dt,
      template = RiboCrypt:::covPanelColumnsPlotlyTemplate()
    )
  )

  expect_true(inherits(p, "plotly"))
  built <- plotly::plotly_build(p)
  expect_equal(unname(built$x$layout$xaxis$range), c(0, 1))
  expect_identical(p$x$source, "mb_top")
})

test_that("annotation_track_allsamples forwards custom regions", {
  display_range <- GenomicRanges::GRangesList(
    tx = GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 10), "+")
  )
  annotation <- GenomicRanges::GRangesList(
    tx = GenomicRanges::GRanges("chr1", IRanges::IRanges(3, 8), "+")
  )
  translons <- GenomicRanges::GRangesList(
    T1 = GenomicRanges::GRanges("chr1", IRanges::IRanges(4, 6), "+")
  )
  captured <- NULL

  testthat::local_mocked_bindings(
    annotation_controller = function(df, display_range, annotation, annotation_names = NULL,
                                     leader_extension, trailer_extension, viewMode) {
      list(display_range = display_range, annotation = annotation)
    },
    createGeneModelPanel = function(display_range, annotation, tx_annotation,
                                    custom_regions, viewMode, collapse_intron_flank,
                                    frame_colors = "R") {
      captured <<- custom_regions
      list(data.table::data.table(), numeric())
    },
    geneModelPanelPlotly = function(dt, template = NULL) plotly::plot_ly(),
    .package = "RiboCrypt"
  )

  p <- RiboCrypt:::annotation_track_allsamples(
    df = NULL,
    id = "tx",
    display_range = display_range,
    annotation = annotation,
    tx_annotation = annotation,
    custom_regions = translons,
    viewMode = "tx",
    collapse_intron_flank = 100
  )

  expect_true(inherits(p, "plotly"))
  expect_identical(captured, translons)
})

test_that("annotation_track_allsamples forwards shared gene model template", {
  display_range <- GenomicRanges::GRangesList(
    tx = GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 10), "+")
  )
  annotation <- GenomicRanges::GRangesList(
    tx = GenomicRanges::GRanges("chr1", IRanges::IRanges(3, 8), "+")
  )
  template <- plotly::plot_ly() %>% plotly::layout(paper_bgcolor = "pink")
  captured_template <- NULL

  testthat::local_mocked_bindings(
    annotation_controller = function(df, display_range, annotation, annotation_names = NULL,
                                     leader_extension, trailer_extension, viewMode) {
      list(display_range = display_range, annotation = annotation)
    },
    createGeneModelPanel = function(display_range, annotation, tx_annotation,
                                    custom_regions, viewMode, collapse_intron_flank,
                                    frame_colors = "R") {
      list(data.table::data.table(), numeric())
    },
    geneModelPanelPlotly = function(dt, template = NULL) {
      captured_template <<- template
      plotly::plot_ly()
    },
    .package = "RiboCrypt"
  )

  p <- RiboCrypt:::annotation_track_allsamples(
    df = NULL,
    id = "tx",
    display_range = display_range,
    annotation = annotation,
    tx_annotation = annotation,
    custom_regions = NULL,
    viewMode = "tx",
    collapse_intron_flank = 100,
    templates = list(gene_model_panel_plotly = template)
  )

  expect_true(inherits(p, "plotly"))
  expect_identical(captured_template, template)
})

test_that("sync_megabrowser_x_shiny resets synced tracks to explicit x range on autorange", {
  calls <- list()
  testthat::local_mocked_bindings(
    plotlyProxy = function(outputId, session) structure(list(id = outputId), class = "plotly_proxy"),
    plotlyProxyInvoke = function(p, method, ...) {
      calls[[length(calls) + 1L]] <<- list(id = p$id, method = method, args = list(...))
      p
    },
    .package = "plotly"
  )

  RiboCrypt:::sync_megabrowser_x_shiny(
    ed = list("xaxis.autorange" = TRUE),
    session = NULL,
    sync_sidebar = FALSE,
    x_reset_range = c(1, 500)
  )

  expect_length(calls, 2)
  expect_equal(vapply(calls, `[[`, character(1), "id"), c("mb_top_summary", "mb_bottom_gene"))
  expect_true(all(vapply(calls, function(call) identical(call$method, "relayout"), logical(1))))
  expect_true(all(vapply(calls, function(call) identical(call$args[[1]][["xaxis.range"]], c(1, 500)), logical(1))))
  expect_true(all(vapply(calls, function(call) identical(call$args[[1]][["xaxis.autorange"]], FALSE), logical(1))))
})

test_that("sync_megabrowser_x_shiny keeps sidebar y zoom aligned without reversal", {
  calls <- list()
  testthat::local_mocked_bindings(
    plotlyProxy = function(outputId, session) structure(list(id = outputId), class = "plotly_proxy"),
    plotlyProxyInvoke = function(p, method, ...) {
      calls[[length(calls) + 1L]] <<- list(id = p$id, method = method, args = list(...))
      p
    },
    .package = "plotly"
  )

  RiboCrypt:::sync_megabrowser_x_shiny(
    ed = list("yaxis.range[0]" = 8.4, "yaxis.range[1]" = 4.6),
    session = NULL,
    sync_tracks = FALSE,
    sync_sidebar = TRUE,
    y_max = 8,
    y_reversed = TRUE
  )

  expect_length(calls, 1)
  expect_identical(calls[[1]]$id, "d")
  expect_identical(calls[[1]]$method, "relayout")
  expect_equal(calls[[1]]$args[[1]]$yaxis$range, c(4.5, 0.5))
  expect_identical(calls[[1]]$args[[1]]$yaxis$autorange, FALSE)
})

test_that("addMegabrowserDoubleClickReset attaches a double-click reset hook", {
  p <- plotly::plot_ly(x = 1:10, y = 1:10, type = "scatter", mode = "lines")
  p <- RiboCrypt:::addMegabrowserDoubleClickReset(
    p,
    reset_range = c(1, 100),
    peer_ids = c("plot-a", "plot-b"),
    reset_layout = list(
      "xaxis.range" = c(1, 100),
      "xaxis.autorange" = FALSE,
      "yaxis.range" = c(10.5, 0.5),
      "yaxis.autorange" = FALSE
    ),
    peer_reset_layout = list(
      "xaxis.range" = c(1, 100),
      "xaxis.autorange" = FALSE
    )
  )

  expect_true(length(p$jsHooks$render) >= 1)
  render_code <- paste(vapply(p$jsHooks$render, `[[`, character(1), "code"), collapse = "\n")
  expect_match(render_code, "plotly_doubleclick")
  expect_match(render_code, "addEventListener\\('dblclick'")
  expect_match(render_code, "reset_layout")
  expect_match(render_code, "peer_reset_layout")
  expect_match(render_code, "peer_ids")
})

test_that("get_meta_browser_plot returns plotly heatmap for plotly type", {
  mat <- matrix(1:20, nrow = 5, ncol = 4)
  table <- data.table::data.table(mat)
  setnames(table, new = paste0("lib", seq(4)))
  data.table::setattr(table, "ratio", 1L)
  km <- stats::kmeans(t(mat), centers = 2)
  data.table::setattr(table, "km", km)
  data.table::setattr(table, "row_order_list", list("1" = seq(4)))


  p <- RiboCrypt:::get_meta_browser_plot(table, color_theme = "default (White-Blue)", plotType = "plotly")
  expect_true(inherits(p, "plotly"))
  expect_identical(p$x$config$doubleClick, FALSE)
  built <- plotly::plotly_build(p)
  expect_equal(unname(built$x$layout$xaxis$range), c(1, 5))
  expect_identical(built$x$layout$xaxis$autorange, FALSE)
  expect_equal(unname(built$x$layout$yaxis$range), c(0.5, 4.5))
  expect_identical(built$x$layout$yaxis$autorange, FALSE)
  expect_identical(built$x$layout$dragmode, "zoom")
  expect_equal(unname(built$x$data[[1]]$z[, 1]), c(1, 6, 11, 16))
})

test_that("get_meta_browser_plot uses original-coordinate x range for binned tables", {
  mat <- matrix(seq_len(600), nrow = 150, ncol = 4)
  table <- data.table::data.table(mat)
  setnames(table, new = paste0("lib", seq(4)))
  data.table::setattr(table, "ratio", 9L)
  data.table::setattr(table, "summary_cov", data.table::data.table(
    count = seq_len(1350),
    position = seq_len(1350),
    frame = factor((seq_len(1350) - 1L) %% 3)
  ))
  km <- stats::kmeans(t(as.matrix(table)), centers = 2)
  data.table::setattr(table, "km", km)
  data.table::setattr(table, "row_order_list", list("1" = seq(2)))

  p <- RiboCrypt:::get_meta_browser_plot(table, color_theme = "default (White-Blue)", plotType = "plotly")
  built <- plotly::plotly_build(p)

  expect_equal(unname(built$x$layout$xaxis$range), c(1, 1350))
  expect_equal(utils::head(built$x$data[[1]]$x, 2), c(1, 10))
})

test_that("get_meta_browser_plot reuses shared heatmap template", {
  mat <- matrix(1:20, nrow = 5, ncol = 4)
  table <- data.table::data.table(mat)
  data.table::setnames(table, new = paste0("lib", seq(4)))
  data.table::setattr(table, "ratio", 1L)
  km <- stats::kmeans(t(mat), centers = 2)
  data.table::setattr(table, "km", km)
  data.table::setattr(table, "row_order_list", list("1" = seq(4)))
  template <- RiboCrypt:::covPanelHeatmapPlotlyTemplate()

  p <- RiboCrypt:::get_meta_browser_plot(
    table,
    color_theme = "default (White-Blue)",
    plotType = "plotly",
    template = RiboCrypt:::megabrowserHeatmapPlotlyTemplate()
  )

  expect_true(inherits(p, "plotly"))
  expect_identical(p$x$data[[1]]$type, "heatmapgl")
  expect_equal(unname(p$x$data[[1]]$x), 1:5)
  expect_equal(dim(p$x$data[[1]]$z), c(4, 5))
})

test_that("get_meta_browser_plot plotly row order matches reversed sidebar order", {
  mat <- matrix(1:20, nrow = 5, ncol = 4)
  table <- data.table::data.table(mat)
  data.table::setnames(table, new = paste0("lib", seq(4)))
  data.table::setattr(table, "ratio", 1L)
  km <- stats::kmeans(t(mat), centers = 2)
  data.table::setattr(table, "km", km)
  data.table::setattr(table, "row_order_list", list("1" = c(2L, 4L, 1L, 3L)))

  p <- RiboCrypt:::get_meta_browser_plot(
    table,
    color_theme = "default (White-Blue)",
    plotType = "plotly",
    template = RiboCrypt:::megabrowserHeatmapPlotlyTemplate()
  )

  expect_equal(unname(p$x$data[[1]]$z[, 1]), c(6, 16, 1, 11))
})

test_that("mb_mid_plot_shiny registers interactive megabrowser source", {
  p <- plotly::plot_ly(
    x = 1:5,
    y = 1:5,
    z = matrix(seq_len(25), nrow = 5),
    type = "heatmapgl"
  ) %>% plotly::layout(yaxis = list(range = c(5.5, 0.5)))

  p <- RiboCrypt:::mb_mid_plot_shiny(p, "plotly")

  expect_identical(p$x$source, "mb_mid")
})

test_that("get_meta_browser_plot preserves matrix heatmap colorscale with shared template", {
  mat <- matrix(1:20, nrow = 5, ncol = 4)
  table <- data.table::data.table(mat)
  data.table::setnames(table, new = paste0("lib", seq(4)))
  data.table::setattr(table, "ratio", 1L)
  km <- stats::kmeans(t(mat), centers = 2)
  data.table::setattr(table, "km", km)
  data.table::setattr(table, "row_order_list", list("1" = seq(4)))
  p <- RiboCrypt:::get_meta_browser_plot(
    table,
    color_theme = "Matrix (black,green,red)",
    plotType = "plotly",
    template = RiboCrypt:::megabrowserHeatmapPlotlyTemplate()
  )

  built <- plotly::plotly_build(p)
  colorscale <- built$x$data[[1]]$colorscale
  expect_identical(colorscale[[1, 2]], "#000000")
  expect_true(any(colorscale[, 2] == "#2CFA1F"))
  expect_true(any(colorscale[, 2] == "#FF2400"))
})

test_that("compute_collection_table_grouping groups metadata with fallback enrichment term", {
  metadata <- data.table::data.table(
    study_accession = c(rep("PRJNA100001", 4), rep("PRJNA100002", 4)),
    Run = c("SRR1001", "SRR1002", "SRR1003", "SRR1004", "SRR2001", "SRR2002", "SRR2003", "SRR2004"),
    BioProject = c(rep("PRJNA100001", 4), rep("PRJNA100002", 4)),
    TISSUE = c("brain", "brain", "brain", "brain", "heart", "heart", "heart", "heart"),
    CELL_LINE = c("CL1", "CL1", "CL2", "CL2", "CL3", "CL3", "CL4", "CL4"),
    CONDITION = c("ctrl", "drug", "ctrl", "drug", "ctrl", "drug", "ctrl", "drug"),
    GENE = c("GENE1", "GENE1", "GENE2", "GENE2", "GENE1", "GENE1", "GENE2", "GENE2")
  )

  table <- data.table::as.data.table(cbind(
    SRR1001 = c(130, 15, 12, 10),   # brain_ctrl
    SRR1002 = c(14, 132, 11, 10),   # brain_drug
    SRR1003 = c(128, 16, 13, 9),    # brain_ctrl
    SRR1004 = c(15, 129, 10, 11),   # brain_drug
    SRR2001 = c(12, 10, 131, 14),   # heart_ctrl
    SRR2002 = c(10, 12, 15, 130),   # heart_drug
    SRR2003 = c(13, 9, 129, 16),    # heart_ctrl
    SRR2004 = c(11, 13, 14, 128)    # heart_drug
  ))
  valid_libs <- rep(TRUE, ncol(table))
  names(valid_libs) <- colnames(table)
  data.table::setattr(table, "valid_libs", valid_libs)

  grouping_args <- list(
    metadata = metadata,
    df = NULL,
    metadata_field = c("TISSUE", "CELL_LINE", "CONDITION", "GENE"),
    table = table,
    ratio_interval = NULL,
    group_on_tx_tpm = NULL,
    decreasing_order = FALSE,
    enrichment_term = "Cluster"
  )

  res <- do.call(RiboCrypt:::compute_collection_table_grouping, grouping_args)

  expect_length(res, 8)
  expect_equal(as.character(res), c("brain", "brain", "brain", "brain", "heart", "heart", "heart", "heart"))
  expect_equal(names(res), c("SRR1001", "SRR1002", "SRR1003", "SRR1004", "SRR2001", "SRR2002", "SRR2003", "SRR2004"))
  expect_equal(attr(res, "xlab"), "TISSUE")
  expect_equal(attr(res, "meta_order"), 1:8)

  other_cols <- attr(res, "other_columns")
  expect_true(all(c("CELL_LINE", "CONDITION", "GENE", "Run") %in% colnames(other_cols)))
  expect_equal(nrow(other_cols), 8)

  run_ids <- attr(res, "runIDs")
  expect_true(all(c("Run", "BioProject") %in% colnames(run_ids)))
  expect_equal(nrow(run_ids), 8)

  # On Gene:
  grouping_args$metadata_field <- c("GENE", "TISSUE", "CELL_LINE", "CONDITION")
  res <- do.call(RiboCrypt:::compute_collection_table_grouping, grouping_args)
  expect_equal(attr(res, "meta_order"), c(1, 2, 5, 6, 3, 4, 7, 8))
  expect_equal(
    names(res),
    c("SRR1001", "SRR1002", "SRR2001", "SRR2002", "SRR1003", "SRR1004", "SRR2003", "SRR2004")
  )
  expect_equal(
    as.character(res),
    c("GENE1", "GENE1", "GENE1", "GENE1", "GENE2", "GENE2", "GENE2", "GENE2")
  )
  expect_equal(attr(res, "xlab"), "GENE")

  ## Test clustering works on the groups
  data.table::setcolorder(table, attr(res, "meta_order"))
  set.seed(1)
  table <- RiboCrypt:::clustering_megabrowser(table, clusters = 4)
  table_with_meta <- list(table = table, metadata_field = res)

  meta_and_clusters <- RiboCrypt:::allsamples_metadata_clustering(
    table_with_meta,
    grouping_args$enrichment_term
  )
  expect_true(is.list(meta_and_clusters))
  expect_true(all(c("meta", "enrich_dt") %in% names(meta_and_clusters)))
  expect_equal(nrow(meta_and_clusters$meta), 8)
  expect_true(all(c("Run", "grouping", "order", "TISSUE", "CELL_LINE", "CONDITION", "cluster", "index") %in%
                    colnames(meta_and_clusters$meta)))
  expect_equal(sort(unique(meta_and_clusters$meta$grouping)), c("GENE1", "GENE2"))
  expect_equal(attr(meta_and_clusters$meta, "xlab"), "GENE Clusters (K-means)")
  expect_equal(attr(meta_and_clusters$meta, "ylab"), "Enrichment")
  expect_true(is.data.frame(meta_and_clusters$enrich_dt))
  expect_equal(rownames(meta_and_clusters$enrich_dt), c("GENE1", "GENE2"))
  expect_equal(meta_and_clusters$meta$cluster, rep(seq(4, 1), each = 2))

  meta_tbl <- RiboCrypt:::allsamples_meta_table(meta_and_clusters)
  expect_equal(nrow(meta_tbl), 8)
  expect_equal(meta_tbl$cluster, as.character(rep(seq(4), each = 2)))
  expect_true(all(c("Run", "BioProject", "grouping", "cluster", "TISSUE", "CELL_LINE", "CONDITION") %in%
                    colnames(meta_tbl)))
  expect_false(any(c("index", "order") %in% colnames(meta_tbl)))
  expect_true(any(meta_tbl$Run == "SRR1001" & meta_tbl$BioProject == "PRJNA100001" & meta_tbl$TISSUE == "brain" &
                    meta_tbl$CELL_LINE == "CL1" & meta_tbl$CONDITION == "ctrl"))

})

test_that("normalize_collection returns a matrix and preserves megabrowser attributes", {
  table <- data.table::data.table(
    SRR1 = c(10, 20, 30),
    SRR2 = c(5, 10, 15)
  )
  valid_libs <- c(SRR1 = TRUE, SRR2 = TRUE)
  data.table::setattr(table, "valid_libs", valid_libs)

  res <- RiboCrypt:::normalize_collection(
    table,
    normalization = "transcriptNormalized",
    lib_sizes = c(SRR1 = 1000, SRR2 = 2000),
    kmer = 1L
  )

  expect_true(is.matrix(res))
  expect_equal(colnames(res), c("SRR1", "SRR2"))
  expect_identical(attr(res, "valid_libs"), valid_libs)
  expect_identical(attr(res, "ratio"), 1L)
  expect_true(is.data.frame(attr(res, "summary_cov")))
  expect_equal(attr(res, "summary_cov")$position, 1:3)
})

test_that("matrix collections still support grouping and clustering", {
  metadata <- data.table::data.table(
    Run = c("SRR1", "SRR2", "SRR3", "SRR4"),
    BioProject = c("P1", "P1", "P2", "P2"),
    TISSUE = c("brain", "brain", "heart", "heart"),
    CONDITION = c("ctrl", "drug", "ctrl", "drug")
  )
  table <- matrix(
    c(
      100, 90, 10, 12,
      95, 85, 15, 14,
      12, 10, 110, 100
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(NULL, metadata$Run)
  )
  attr(table, "valid_libs") <- setNames(rep(TRUE, ncol(table)), colnames(table))

  grouping <- RiboCrypt:::compute_collection_table_grouping(
    metadata = metadata,
    df = NULL,
    metadata_field = c("TISSUE", "CONDITION"),
    table = table,
    ratio_interval = NULL,
    group_on_tx_tpm = NULL,
    decreasing_order = FALSE,
    enrichment_term = "TISSUE"
  )

  ordered_table <- table[, attr(grouping, "meta_order"), drop = FALSE]
  ordered_table <- RiboCrypt:::clustering_megabrowser(ordered_table, clusters = 2)

  expect_true(is.matrix(ordered_table))
  expect_length(attr(ordered_table, "row_order_list"), 2)
  expect_identical(attr(grouping, "xlab"), "TISSUE")
})

test_that("matrix column subsetting preserves megabrowser attributes", {
  table <- matrix(
    c(1, 2, 3, 4, 5, 6),
    nrow = 3,
    dimnames = list(NULL, c("SRR1", "SRR2"))
  )
  attr(table, "summary_cov") <- data.table::data.table(
    count = c(3, 7, 11),
    position = 1:3,
    frame = factor(c(0, 1, 2))
  )
  attr(table, "ratio") <- 1L
  attr(table, "valid_libs") <- c(SRR1 = TRUE, SRR2 = TRUE)

  res <- RiboCrypt:::subset_collection_columns(table, c(2, 1))

  expect_true(is.matrix(res))
  expect_equal(colnames(res), c("SRR2", "SRR1"))
  expect_true(is.data.frame(attr(res, "summary_cov")))
  expect_equal(attr(res, "summary_cov")$position, 1:3)
  expect_identical(attr(res, "ratio"), 1L)
  expect_identical(attr(res, "valid_libs"), c(SRR1 = TRUE, SRR2 = TRUE))
})

test_that("allsamples_sidebar_plotly reuses shared numeric and categorical templates", {
  meta <- data.table::data.table(
    Run = c("SRR1", "SRR2", "SRR3"),
    grouping = c("A", "B", "C"),
    order = 1:3,
    index = 1:3,
    score = c(10, 20, 30),
    tissue = c("brain", "heart", "brain")
  )
  attr(meta, "xlab") <- "grouping"

  p <- RiboCrypt:::allsamples_sidebar_plotly(
    meta,
    templates = list(
      allsamples_sidebar_numeric_plotly = RiboCrypt:::allsamplesSidebarNumericPlotlyTemplate(),
      allsamples_sidebar_categorical_plotly = RiboCrypt:::allsamplesSidebarCategoricalPlotlyTemplate()
    )
  )

  expect_true(inherits(p, "plotly"))
  built <- plotly::plotly_build(p)
  expect_identical(built$x$data[[1]]$type, "scatter")
  expect_equal(as.numeric(built$x$data[[1]]$x), c(10, 20, 30))
  expect_equal(as.integer(built$x$data[[1]]$y), 1:3)
  expect_identical(built$x$data[[2]]$type, "heatmap")
  expect_equal(as.integer(built$x$data[[2]]$y), 1:3)
  expect_equal(dim(built$x$data[[2]]$z), c(3, 1))
  expect_gt(length(unique(built$x$data[[2]]$colorscale[, 2])), 1)
})

test_that("allsamples_sidebar_plotly reuses shared cluster label template", {
  meta <- data.table::data.table(
    Run = c("SRR1", "SRR2", "SRR3", "SRR4"),
    grouping = c("A", "A", "B", "B"),
    order = 1:4,
    index = 1:4,
    cluster = c(1, 1, 2, 2),
    tissue = c("brain", "brain", "heart", "heart")
  )
  attr(meta, "xlab") <- "grouping"

  p <- RiboCrypt:::allsamples_sidebar_plotly(
    meta,
    templates = list(
      allsamples_sidebar_categorical_plotly = RiboCrypt:::allsamplesSidebarCategoricalPlotlyTemplate(),
      allsamples_sidebar_cluster_label_plotly = RiboCrypt:::allsamplesSidebarClusterLabelPlotlyTemplate()
    )
  )

  expect_true(inherits(p, "plotly"))
  built <- plotly::plotly_build(p)
  expect_identical(built$x$data[[1]]$type, "scatter")
  expect_identical(built$x$data[[2]]$type, "heatmap")
  expect_length(built$x$layout$annotations, 2)
  expect_equal(vapply(built$x$layout$annotations, `[[`, character(1), "text"), c("1", "2"))
})
