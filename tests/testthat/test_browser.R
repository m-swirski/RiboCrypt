df <- ORFik::ORFik.template.experiment()[9:10, ]
tx <- ORFik::loadRegion(df, "mrna")
cds <- ORFik::loadRegion(df, "cds")

tx_unlisted <- unlist(tx[1])
cds_unlisted <- unlist(cds[1])
# Overlapping uORF
T1 <- GenomicRanges::GRanges(
  seqnames = GenomicRanges::seqnames(tx_unlisted)[1],
  ranges = IRanges::IRanges(start = c(start(cds_unlisted) - 16, start(cds_unlisted)),
                            end = c(end(tx_unlisted[1]), start(cds_unlisted) + 13)),
  strand = BiocGenerics::strand(tx_unlisted)[1]
)
# Non overlapping uORF
T2 <- GenomicRanges::GRanges(
  seqnames = GenomicRanges::seqnames(tx_unlisted)[1],
  ranges = IRanges::IRanges(start = c(start(tx_unlisted[1]) + 30),
                            end = c(start(tx_unlisted[1]) + 60)),
  strand = BiocGenerics::strand(tx_unlisted)[1]
)

# Helper functions
make_bottom_panel_test_controls <- function(viewMode = FALSE,
                                            df, tx, cds,
                                            tx_id = names(tx)[1],
                                            collapsed_introns_width = 0,
                                            genomic_string = NULL,
                                            extendLeaders = 0,
                                            extendTrailers = 0,
                                            frame_colors = "R") {

  display_region <- RiboCrypt:::observed_tx_annotation(tx_id, function() tx)
  annotation <- RiboCrypt:::observed_cds_annotation_internal(tx_id, cds, TRUE)
  tx_annotation <- RiboCrypt:::observed_cds_annotation_internal(tx_id, tx, TRUE)

  if (isTRUE(viewMode) && collapsed_introns_width > 0) {
    tx_annotation <- tx_annotation[tx_annotation %over% ORFik::flankPerGroup(display_region)]
    display_region_gr <- GenomicRanges::reduce(ORFik::unlistGrl(tx_annotation))
    display_region <- ORFik::groupGRangesBy(
      display_region_gr,
      rep(names(display_region), length(display_region_gr))
    )
  }
  display_region <- genomic_string_to_grl(genomic_string, display_region,
                                          max_size = 1e6, viewMode,
                                          0,
                                          0,
                                          collapsed_introns_width)

  controls <- list(
    dff = df,
    display_region = display_region,
    annotation = annotation,
    extendLeaders = extendLeaders,
    extendTrailers = extendTrailers,
    viewMode = viewMode,
    custom_sequence = "",
    customRegions = GenomicRanges::GRangesList(),
    tx_annotation = tx_annotation,
    collapsed_introns_width = collapsed_introns_width,
    frame_colors = frame_colors,
    gg_theme = RiboCrypt:::gg_theme_template(),
    phyloP = FALSE,
    mapability = FALSE,
    is_cellphone = FALSE
  )
  list(
    controls = function() controls,
    session = list(ns = function(x) x)
  )
}

make_browser_track_test_fixture <- function(frames_type, kmers, df, tx, cds) {
  df <- df[1,] # 1 is 9
  tx_id <- names(tx)[1]
  reads <- filepath(df, "bigwig")

  controls <- list(
    dff = df,
    display_region = RiboCrypt:::observed_tx_annotation(tx_id, function() tx),
    annotation = RiboCrypt:::observed_cds_annotation_internal(tx_id, cds, TRUE),
    extendLeaders = 0,
    extendTrailers = 0,
    viewMode = FALSE,
    custom_sequence = "",
    customRegions = GenomicRanges::GRangesList(),
    tx_annotation = RiboCrypt:::observed_cds_annotation_internal(tx_id, tx, TRUE),
    collapsed_introns_width = 0,
    frame_colors = "R",
    gg_theme = RiboCrypt:::gg_theme_template(),
    phyloP = FALSE,
    mapability = FALSE,
    is_cellphone = FALSE,
    reads = reads,
    withFrames = TRUE,
    frames_type = frames_type,
    colors = NULL,
    kmerLength = kmers,
    log_scale = FALSE,
    summary_track = FALSE,
    summary_track_type = frames_type,
    export_format = "svg",
    zoom_range = numeric(0),
    frames_subset = "all"
  )

  list(
    controls = function() controls,
    session = list(ns = function(x) x)
  )
}

make_collapsed_overlap_fixture <- function() {
  display_tx <- GenomicRanges::GRangesList(
    txA = GenomicRanges::GRanges(
      "chr1",
      IRanges::IRanges(c(100, 500, 900), c(199, 599, 999)),
      "+"
    )
  )

  cds_overlap <- GenomicRanges::GRangesList(
    txA = GenomicRanges::GRanges(
      "chr1",
      IRanges::IRanges(c(100, 500, 900), c(199, 599, 999)),
      "+"
    ),
    txB = GenomicRanges::GRanges(
      "chr1",
      IRanges::IRanges(c(450, 900), c(520, 950)),
      "+"
    )
  )

  tx_overlap <- GenomicRanges::GRangesList(
    txA = display_tx[[1]],
    txB = GenomicRanges::GRanges(
      "chr1",
      IRanges::IRanges(c(450, 520, 900, 980), c(520, 560, 950, 999)),
      "+"
    )
  )

  list(
    display_range = ORFik::exonsWithPseudoIntronsPerGroup(display_tx, 30),
    annotation = cds_overlap,
    tx_annotation = tx_overlap
  )
}


test_that("input_to_list drops ignored inputs and adds user info", {
  input <- shiny::reactiveValues(
    gene = "GENE1",
    tx = "TX1",
    go = 1,
    toggle_settings = 0,
    select_all_btn = 0,
    c__shinyjquiBookmarkState__resizable = NULL,
    c_is_resizing = FALSE,
    c_size = NULL,
    browser_plot__shinyjquiBookmarkState__resizable = NULL,
    browser_plot_is_resizing = FALSE,
    browser_plot_size = NULL
  )
  res <- shiny::isolate(RiboCrypt:::input_to_list(
    input,
    list(id = "user-1", browser = "firefox")
  ))
  expect_true("gene" %in% names(res))
  expect_true("tx" %in% names(res))
  expect_false("go" %in% names(res))
  expect_false("toggle_settings" %in% names(res))
  expect_false("select_all_btn" %in% names(res))
  expect_false("c__shinyjquiBookmarkState__resizable" %in% names(res))
  expect_false("c_is_resizing" %in% names(res))
  expect_false("c_size" %in% names(res))
  expect_false("browser_plot__shinyjquiBookmarkState__resizable" %in% names(res))
  expect_false("browser_plot_is_resizing" %in% names(res))
  expect_false("browser_plot_size" %in% names(res))
  expect_true("user.info.browser" %in% names(res))
  expect_equal(res[["user.info.browser"]], "firefox")
})

test_that("go_when_input_is_ready triggers kickoff when inputs match", {
  browser_options <- c(
    plot_on_start = "TRUE",
    default_gene = "GENE1",
    default_isoform = "TX1",
    default_libs = "A,B"
  )
  input <- list(
    gene = "GENE1",
    tx = "TX1",
    library = c("A", "B")
  )
  libs <- shiny::reactiveVal(c("A", "B"))
  fired <- shiny::reactiveVal(FALSE)
  kickoff <- shiny::reactiveVal(FALSE)

  shiny::isolate(RiboCrypt:::go_when_input_is_ready(input, browser_options, fired, kickoff, libs))
  expect_true(shiny::isolate(isTRUE(fired())))
  expect_true(shiny::isolate(isTRUE(kickoff())))
})

test_that("go_when_input_is_ready does not trigger when plot_on_start is FALSE", {
  browser_options <- c(
    plot_on_start = "FALSE",
    default_gene = "GENE1",
    default_isoform = "TX1",
    default_libs = "A,B"
  )
  input <- list(
    gene = "GENE1",
    tx = "TX1",
    library = c("A", "B")
  )
  libs <- shiny::reactiveVal(c("A", "B"))
  fired <- shiny::reactiveVal(FALSE)
  kickoff <- shiny::reactiveVal(FALSE)

  shiny::isolate(RiboCrypt:::go_when_input_is_ready(input, browser_options, fired, kickoff, libs))
  expect_true(shiny::isolate(isTRUE(fired())))
  expect_false(shiny::isolate(isTRUE(kickoff())))
})

test_that("browser_ui returns a shiny tag object", {
  all_exp <- data.frame(
    organism = c("org1", "org2"),
    name = c("exp1", "exp2"),
    stringsAsFactors = FALSE
  )
  browser_options <- c(
    default_gene = "GENE1",
    default_isoform = "TX1",
    default_libs = "A|B",
    default_view_mode = "tx",
    collapsed_introns_width = "100",
    full_annotation = "FALSE",
    translons = "FALSE",
    translons_transcode = "FALSE",
    hide_settings = "TRUE",
    default_frame_type = "all",
    default_kmer = "1"
  )
  gene_names_init <- data.frame(label = "GENE1", stringsAsFactors = FALSE)
  libs <- c("A", "B")

  ui <- shiny::isolate(RiboCrypt:::browser_ui("test", all_exp, browser_options, gene_names_init, libs))
  expect_true(inherits(ui, "shiny.tag") || inherits(ui, "shiny.tag.list"))
})


test_that("bottom_panel_shiny handles transcript view", {
  controls <- make_bottom_panel_test_controls(viewMode = FALSE, df, tx, cds)

  bottom_panel <- RiboCrypt:::bottom_panel_shiny(controls$controls)

  expect_equal(ORFik::widthPerGroup(bottom_panel$display_range, FALSE), 551)
  expect_equal(unname(bottom_panel$lines), c(141, 446))
  expect_equal(length(bottom_panel$bottom_plots), 3)
  expect_equal(bottom_panel$annotation_layers, 1)
})

test_that("bottom_panel_shiny applies 5' and 3' extensions", {
  original_size <- widthPerGroup(tx[1], FALSE)
  extendLeaders <- 10
  extendTrailers <- 20
  new_size <- original_size + extendLeaders + extendTrailers
  controls <- make_bottom_panel_test_controls(
    viewMode = FALSE, df, tx, cds,
    extendLeaders = extendLeaders,
    extendTrailers = extendTrailers
  )
  bottom_panel <- RiboCrypt:::bottom_panel_shiny(controls$controls)
  expect_equal(ORFik::widthPerGroup(bottom_panel$display_range, FALSE), new_size)


})

test_that("bottom_panel_shiny supports Color_blind frame theme", {
  controls <- make_bottom_panel_test_controls(
    viewMode = FALSE, df, tx, cds,
    frame_colors = "Color_blind"
  )

  bottom_panel <- RiboCrypt:::bottom_panel_shiny(controls$controls)
  cds_cols <- unique(bottom_panel$gene_model_panel_dt[type == "cds"]$cols)

  expect_true(all(cds_cols %in% RiboCrypt:::frame_color_themes("Color_blind", TRUE)))
})

test_that("bottom_panel_shiny handles transcript view with custom region overlapping CDS", {
  controls <- make_bottom_panel_test_controls(viewMode = FALSE, df, tx, cds)
  controls_with_custom <- controls$controls()
  tx_unlisted <- unlist(controls_with_custom$display_region[1])
  cds_unlisted <- unlist(controls_with_custom$annotation[1])
  controls_with_custom$customRegions <- GenomicRanges::GRangesList(
    T1 = T1
  )

  bottom_panel <- RiboCrypt:::bottom_panel_shiny(function() controls_with_custom)

  expect_equal(ORFik::widthPerGroup(bottom_panel$display_range, FALSE), 551)
  expect_equal(length(bottom_panel$bottom_plots), 3)
  expect_equal(bottom_panel$annotation_layers, 2)
})

test_that("bottom_panel_shiny handles transcript view with custom region overlapping CDS
          and one not overlapping any", {
  controls <- make_bottom_panel_test_controls(viewMode = FALSE, df, tx, cds)
  controls_with_custom <- controls$controls()

  controls_with_custom$customRegions <- GenomicRanges::GRangesList(
    T1 = T1,
    T2 = T2
  )

  bottom_panel <- RiboCrypt:::bottom_panel_shiny(function() controls_with_custom)

  expect_equal(ORFik::widthPerGroup(bottom_panel$display_range, FALSE), 551)
  expect_equal(length(bottom_panel$bottom_plots), 3)
  expect_equal(bottom_panel$annotation_layers, 2)
})

test_that("bottom_panel_shiny handles genomic view (non-spliced)", {
  controls <- make_bottom_panel_test_controls(viewMode = TRUE, df, tx, cds)

  bottom_panel <- RiboCrypt:::bottom_panel_shiny(controls$controls)

  expect_equal(ORFik::widthPerGroup(bottom_panel$display_range, FALSE), 551)
  expect_length(unlist(bottom_panel$display_range), 1)
  expect_equal(as.numeric(width(unlist(bottom_panel$display_range))), c(551))
  expect_equal(length(bottom_panel$bottom_plots), 3)
  expect_equal(bottom_panel$annotation_layers, 1)
})

test_that("bottom_panel_shiny handles genomic view (non-spliced)
          with custom region overlapping CDS", {
  controls <- make_bottom_panel_test_controls(viewMode = TRUE, df, tx, cds)
  controls_with_custom <- controls$controls()

  controls_with_custom$customRegions <- GenomicRanges::GRangesList(
    T1 = T1
  )
  bottom_panel <- RiboCrypt:::bottom_panel_shiny(function() controls_with_custom)
  dt <- bottom_panel$gene_model_panel_dt # Layer box
  expect_equal(ORFik::widthPerGroup(bottom_panel$display_range, FALSE), 551)
  expect_length(unlist(bottom_panel$display_range), 1)
  expect_equal(as.numeric(width(unlist(bottom_panel$display_range))), c(551))
  expect_equal(length(bottom_panel$bottom_plots), 3)
  expect_equal(bottom_panel$annotation_layers, 2)
  expect_equal(unique(dt[gene_names == "ENSTTEST10001"]$layers), 1)
  expect_equal(unique(dt[gene_names == "T1"]$layers), 2)
})

test_that("bottom_panel_shiny handles genomic view (spliced) with collapsed introns", {
  tx_spliced <- tx[1]
  end(tx_spliced[1])[[1]][1] <- 350
  controls <- make_bottom_panel_test_controls(
    viewMode = TRUE, df, tx_spliced, cds,
    collapsed_introns_width = 30
  )

  bottom_panel <- RiboCrypt:::bottom_panel_shiny(controls$controls)

  expect_equal(ORFik::widthPerGroup(bottom_panel$display_range, FALSE), 516)
  expect_length(unlist(bottom_panel$display_range), 2)
  expect_equal(as.numeric(width(unlist(bottom_panel$display_range))), c(75, 441))
  expect_equal(length(bottom_panel$bottom_plots), 3)
  expect_equal(bottom_panel$annotation_layers, 1)
})

test_that("bottom_panel_shiny keeps overlapping CDS exons in collapsed genomic view", {
  fixture <- make_collapsed_overlap_fixture()
  controls <- list(
    dff = df,
    display_region = fixture$display_range,
    annotation = fixture$annotation,
    extendLeaders = 0,
    extendTrailers = 0,
    viewMode = TRUE,
    custom_sequence = "",
    customRegions = GenomicRanges::GRangesList(),
    tx_annotation = fixture$tx_annotation,
    collapsed_introns_width = 30,
    frame_colors = "R",
    gg_theme = RiboCrypt:::gg_theme_template(),
    phyloP = FALSE,
    mapability = FALSE,
    is_cellphone = FALSE
  )

  bottom_panel <- suppressWarnings(RiboCrypt:::bottom_panel_shiny(function() controls))
  main_cds <- bottom_panel$gene_model_panel_dt[gene_names == "txA" & type == "cds"]
  overlapping_cds <- bottom_panel$gene_model_panel_dt[gene_names == "txB" & type == "cds"]

  # TODO: Decide if this is smart or not that RIboCrypt just removes txB here.
  expect_true("rect_starts" %in% colnames(overlapping_cds))
  expect_equal(nrow(overlapping_cds), 0)
  expect_equal(main_cds$rect_starts, c(1, 161, 321))
  expect_equal(main_cds$rect_ends, c(100,260,420))
})

test_that("bottom_panel_shiny handles genomic string input on non tx region", {
  genomic_string <- "chr1:1-300:+"
  # Tx is ignored when genomic string is given
  controls <- make_bottom_panel_test_controls(
    viewMode = TRUE, df, tx, cds,
    genomic_string = genomic_string
  )

  bottom_panel <- RiboCrypt:::bottom_panel_shiny(controls$controls)

  expect_equal(ORFik::widthPerGroup(bottom_panel$display_range, FALSE), 300)
  expect_length(unlist(bottom_panel$display_range), 1)
  expect_equal(as.numeric(width(unlist(bottom_panel$display_range))), c(300))
  expect_equal(length(bottom_panel$bottom_plots), 3)
  expect_equal(bottom_panel$annotation_layers, 1)
})


test_that("browser_track_panel_shiny handles area tracks with 9-mers", {
  fixture <- make_browser_track_test_fixture(frames_type = "area", kmers = 9,
                                             df, tx, cds)

  bottom_panel <- RiboCrypt:::bottom_panel_shiny(fixture$controls)
  plot <- RiboCrypt:::browser_track_panel_shiny(
    fixture$controls,
    bottom_panel,
    fixture$session
  )

  expect_s3_class(plot, "plotly")
  expect_true(inherits(plot, "htmlwidget"))
  expect_gt(length(plot$x$data), 0)
})

test_that("browser_track_panel_shiny handles column tracks with single-nucleotide bins", {
  fixture <- make_browser_track_test_fixture(frames_type = "columns", kmers = 1,
                                             df, tx, cds)

  bottom_panel <- RiboCrypt:::bottom_panel_shiny(fixture$controls)
  plot <- RiboCrypt:::browser_track_panel_shiny(
    fixture$controls,
    bottom_panel,
    fixture$session
  )

  expect_s3_class(plot, "plotly")
  expect_true(inherits(plot, "htmlwidget"))
  expect_gt(length(plot$x$data), 0)
})

test_that("browser_track_panel_shiny handles overlapping ranges", {
  fixture <- make_browser_track_test_fixture(frames_type = "columns", kmers = 1,
                                             df, tx, cds)
  controls <- make_bottom_panel_test_controls(viewMode = FALSE, df, tx, cds)
  controls_with_custom <- controls$controls()

  controls_with_custom$customRegions <- GenomicRanges::GRangesList(
    T1 = T1,
    T2 = T2
  )

  bottom_panel <- RiboCrypt:::bottom_panel_shiny(function() controls_with_custom)
  plot <- RiboCrypt:::browser_track_panel_shiny(
    fixture$controls,
    bottom_panel,
    fixture$session
  )

  expect_s3_class(plot, "plotly")
  expect_true(inherits(plot, "htmlwidget"))
  expect_gt(length(plot$x$data), 0)
})

test_that("browser_track_panel_shiny applies zoom_range to the shared x axis", {
  fixture <- make_browser_track_test_fixture(frames_type = "columns", kmers = 1,
                                             df, tx, cds)
  controls_with_zoom <- fixture$controls()
  controls_with_zoom$zoom_range <- c(20, 60)

  bottom_panel <- RiboCrypt:::bottom_panel_shiny(function() controls_with_zoom)
  plot <- RiboCrypt:::browser_track_panel_shiny(
    function() controls_with_zoom,
    bottom_panel,
    fixture$session
  )

  expect_equal(unname(plot$x$layout$xaxis$range), c(20, 60))
})

test_that("browser_track_panel_shiny keeps zoom_range when a custom bottom track is present", {
  fixture <- make_browser_track_test_fixture(frames_type = "columns", kmers = 1,
                                             df, tx, cds)
  controls_with_zoom <- fixture$controls()
  controls_with_zoom$zoom_range <- c(20, 60)

  bottom_panel <- RiboCrypt:::bottom_panel_shiny(function() controls_with_zoom)
  custom_track <- plotly::plot_ly(
    x = 1:100,
    y = rep(1, 100),
    type = "bar",
    showlegend = FALSE
  ) %>%
    plotly::layout(
      xaxis = list(showticklabels = FALSE),
      yaxis = list(title = list(text = "P"), showticklabels = FALSE)
    )
  bottom_panel$custom_bigwig_panels <- list(custom_track)
  bottom_panel$bottom_plots <- RiboCrypt:::bottom_plots_to_plotly(bottom_panel)
  bottom_panel$ncustom <- 1L

  plot <- RiboCrypt:::browser_track_panel_shiny(
    function() controls_with_zoom,
    bottom_panel,
    fixture$session
  )

  xaxis_names <- names(plot$x$layout)[grepl("^xaxis[0-9]*$", names(plot$x$layout))]
  ranges <- lapply(xaxis_names, function(axis_name) plot$x$layout[[axis_name]]$range)

  expect_true(all(vapply(ranges, function(r) identical(unname(r), c(20, 60)), logical(1))))
})

test_that("custom_seq_track_panel_bigwig builds a native plotly bar track", {
  grl <- GenomicRanges::GRangesList(
    tx = GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 4), "+")
  )

  testthat::local_mocked_bindings(
    coveragePerTiling = function(grl, bigwig_path) {
      data.frame(position = 1:4, count = c(0, 2, 1, 3))
    },
    .package = "RiboCrypt"
  )
  testthat::local_mocked_bindings(
    BigWigFile = function(path) grl,
    .package = "rtracklayer"
  )

  plot <- RiboCrypt:::custom_seq_track_panel_bigwig(grl, "dummy.bw", "P")
  built <- plotly::plotly_build(plot)

  expect_s3_class(plot, "plotly")
  expect_true(inherits(plot, "htmlwidget"))
  expect_identical(built$x$data[[1]]$type, "bar")
  expect_identical(built$x$layout$yaxis$title$text, "P")
  expect_false(isTRUE(built$x$data[[1]]$showlegend))
  expect_null(built$x$layout$width)
})

test_that("automateTicksCustomTrack supports native plotly inputs", {
  plot <- plotly::plot_ly(
    x = 1:4,
    y = c(0, 2, 1, 3),
    type = "bar",
    showlegend = FALSE
  )

  styled <- RiboCrypt:::automateTicksCustomTrack(plot)
  built <- plotly::plotly_build(styled)

  expect_s3_class(styled, "plotly")
  expect_true(isTRUE(built$x$layout$yaxis$fixedrange))
  expect_false(isTRUE(built$x$data[[1]]$showlegend))
})

test_that("lineDeSimplify only updates pure line traces", {
  p <- plotly::plot_ly() %>%
    plotly::add_lines(x = 1:3, y = c(1, 2, 3), line = list(color = "red")) %>%
    plotly::add_trace(
      x = c(1, 1, NA, 2, 2),
      y = c(0, 1, NA, 0, 1),
      type = "scatter",
      mode = "lines+markers",
      line = list(color = "white", width = 2),
      marker = list(opacity = 0),
      showlegend = FALSE
    )

  built <- RiboCrypt:::lineDeSimplify(p)

  expect_false(is.null(built$x$data[[1]]$line$simplify))
  expect_identical(built$x$data[[1]]$line$simplify, FALSE)
  expect_null(built$x$data[[2]]$line$simplify)
})

test_that("plotAASeqPanelPlotly hides frame tick labels", {
  hits <- data.table::data.table(
    col = "white",
    pos = c(3L, 9L),
    frames = c(0L, 1L)
  )

  p <- RiboCrypt:::plotAASeqPanelPlotly(hits, Biostrings::DNAString("ATGATGATGATG"))

  expect_false(isTRUE(p$x$layout$yaxis$showticklabels))
})

test_that("browser_plot_final_layout_polish keeps x ticks only on bottom shared axis", {
  p <- plotly::subplot(
    plotly::plot_ly(x = 1:3, y = 1:3, type = "scatter", mode = "lines"),
    plotly::plot_ly(x = 1:3, y = 3:1, type = "scatter", mode = "lines"),
    nrows = 2,
    shareX = TRUE,
    titleX = TRUE
  )

  polished <- RiboCrypt:::browser_plot_final_layout_polish(
    p,
    plot_name = "default",
    display_range = GenomicRanges::GRangesList(
      tx = GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 10), "+")
    ),
    width = NULL,
    height = NULL,
    export.format = "svg",
    plot_title = NULL,
    zoom_range = NULL,
    proportions = c(0.5, 0.5)
  )

  expect_true(isTRUE(polished$x$layout$xaxis$visible))
  expect_true(isTRUE(polished$x$layout$xaxis$showticklabels))
  expect_identical(polished$x$layout$xaxis$title$text, "position [nt]")
  expect_identical(polished$x$layout$legend$y, 0.93)
  expect_identical(polished$x$layout$legend$yanchor, "top")
})

test_that("browser_plot_final_layout_polish applies zoom_range to non-default shared axes", {
  testthat::local_mocked_bindings(
    lineDeSimplify = function(plot) plot,
    .package = "RiboCrypt"
  )

  p <- plotly::plot_ly(x = 1:100, y = 1:100, type = "scatter", mode = "lines")
  p$x$layout$xaxis <- list(anchor = "y", domain = c(0, 1))
  p$x$layout$xaxis2 <- list(anchor = "y2", domain = c(0, 1))
  p$x$layout$yaxis <- list(domain = c(0.4, 1))
  p$x$layout$yaxis2 <- list(domain = c(0, 0.3))

  polished <- RiboCrypt:::browser_plot_final_layout_polish(
    p,
    plot_name = "default",
    display_range = GenomicRanges::GRangesList(
      tx = GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 100), "+")
    ),
    width = NULL,
    height = NULL,
    export.format = "svg",
    plot_title = NULL,
    zoom_range = c(20, 60),
    proportions = c(0.7, 0.3)
  )

  expect_equal(unname(polished$x$layout$xaxis$range), c(20, 60))
  expect_equal(unname(polished$x$layout$xaxis2$range), c(20, 60))
  expect_false(isTRUE(polished$x$layout$xaxis$visible))
  expect_true(isTRUE(polished$x$layout$xaxis2$visible))
  expect_identical(polished$x$layout$xaxis2$title$text, "position [nt]")
})

test_that("lineDeSimplify keeps one legend item per shared frame", {
  profile <- data.table::data.table(
    position = 1:6,
    count = c(1, 2, 3, 2, 1, 0),
    frame = factor(c(0, 1, 2, 0, 1, 2))
  )

  p1 <- RiboCrypt:::createSinglePlot(
    profile, TRUE, "R", NULL, "a", "a", numeric(),
    type = "lines", lib_index = 1, total_libs = 2
  )
  p2 <- RiboCrypt:::createSinglePlot(
    profile, TRUE, "R", NULL, "b", "b", numeric(),
    type = "lines", lib_index = 2, total_libs = 2
  )

  polished <- RiboCrypt:::lineDeSimplify(
    plotly::subplot(list(p1, p2), nrows = 2, shareX = TRUE, titleY = TRUE, titleX = TRUE)
  )

  legend_traces <- Filter(function(tr) isTRUE(tr$showlegend), polished$x$data)
  legend_names <- vapply(legend_traces, function(tr) as.character(if (is.null(tr$name)) "" else tr$name), character(1))

  expect_equal(sort(legend_names[nzchar(legend_names)]), c("0", "1", "2"))
})
