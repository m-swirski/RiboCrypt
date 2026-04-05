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

test_that("browser selection helpers prefer valid species-specific defaults", {
  gene_name_list <- data.table::data.table(
    value = c("TXA1", "TXA2", "TXB1"),
    label = c("GENEA", "GENEA", "GENEB")
  )

  expect_equal(
    RiboCrypt:::resolve_gene_selection(gene_name_list, preferred = "GENEB", fallback = "GENEA"),
    "GENEB"
  )
  expect_equal(
    RiboCrypt:::resolve_gene_selection(gene_name_list, preferred = "MISSING", fallback = "GENEA"),
    "GENEA"
  )
  expect_equal(
    RiboCrypt:::resolve_tx_selection(gene_name_list, gene = "GENEA", preferred = "TXA2", fallback = "TXA1"),
    "TXA2"
  )
  expect_equal(
    RiboCrypt:::resolve_tx_selection(gene_name_list, gene = "GENEA", preferred = "MISSING", fallback = "TXA1"),
    "TXA1"
  )
  expect_equal(
    RiboCrypt:::resolve_tx_selection(gene_name_list, gene = "MISSING", preferred = "TXA1"),
    character()
  )
})

test_that("multiOmicsControllerView supports default selected_libraries", {
  env <- rlang::env(
    reads = list(1, 2),
    withFrames = logical(),
    colors = NULL,
    kmers = NULL,
    ylabels = NULL,
    lib_proportions = NULL,
    annotation_proportions = NULL,
    frames_type = "lines",
    display_sequence = "none",
    viewMode = "tx",
    bottom_panel = list(annotation_layers = 1, ncustom = 0),
    lib_to_annotation_proportions = c(0.8, 0.2),
    summary_track = FALSE
  )

  expect_no_error(evalq(RiboCrypt:::multiOmicsControllerView(), env))
  expect_equal(env$withFrames, c(FALSE, FALSE))
  expect_equal(env$colors, c(1L, 2L))
  expect_equal(env$kmers, c(1, 1))
  expect_equal(env$ylabels, c("1", "2"))
  expect_equal(env$ylabels_full_name, c("1", "2"))
  expect_equal(env$lib_proportions, c(0.5, 0.5))
  expect_equal(env$annotation_proportions, c(0.35, 0.65))
  expect_equal(env$proportions, c(0.4, 0.4, 0.07, 0.13))
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
  expect_s3_class(bottom_panel$seq_nt_panel, "plotly")
  expect_equal(bottom_panel$seq_nt_panel$x$attrs[[1]]$x, c(1, 551))
  expect_equal(bottom_panel$seq_nt_panel$x$attrs[[1]]$y, c(0, 1))
  expect_equal(bottom_panel$seq_nt_panel$x$layoutAttrs[[1]]$xaxis$range, c(1, 551))
  expect_equal(bottom_panel$seq_nt_panel$x$layoutAttrs[[1]]$yaxis$range, c(0, 1))
  expect_false(isTRUE(bottom_panel$seq_nt_panel$x$layoutAttrs[[1]]$xaxis$showticklabels))
  expect_false(isTRUE(bottom_panel$seq_nt_panel$x$layoutAttrs[[1]]$yaxis$showticklabels))
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
  aa_shape_cols <- vapply(bottom_panel$seq_panel$x$layout$shapes, `[[`, character(1), "fillcolor")

  expect_true(all(cds_cols %in% RiboCrypt:::frame_color_themes("Color_blind", FALSE)))
  expect_equal(unname(aa_shape_cols), RiboCrypt:::frame_color_themes("Color_blind", FALSE))
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

  expect_no_warning(
    bottom_panel <- RiboCrypt:::bottom_panel_shiny(controls$controls)
  )

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

test_that("plotAASeqPanelPlotly keeps codon marker lines unsimplified", {
  hits <- data.table::data.table(
    col = c("white", "black"),
    pos = c(3L, 9L),
    frames = c(0L, 1L)
  )

  built <- plotly::plotly_build(
    RiboCrypt:::plotAASeqPanelPlotly(hits, Biostrings::DNAString("ATGATGATGATG"))
  )
  marker_traces <- built$x$data[vapply(built$x$data, function(tr) identical(tr$mode, "lines"), logical(1))]

  expect_length(marker_traces, 2)
  expect_true(all(vapply(marker_traces, function(tr) identical(tr$type, "scatter"), logical(1))))
  expect_true(all(vapply(marker_traces, function(tr) identical(tr$line$simplify, FALSE), logical(1))))
  expect_true(all(vapply(marker_traces, function(tr) identical(tr$line$width, 2), logical(1))))
})

test_that("plotAASeqPanelPlotly reuses template shapes without mutating template", {
  hits <- data.table::data.table(
    col = "white",
    pos = c(3L, 9L),
    frames = c(0L, 1L)
  )
  template <- RiboCrypt:::aaSeqPanelPlotlyTemplate()

  p <- RiboCrypt:::plotAASeqPanelPlotly(
    hits,
    Biostrings::DNAString("ATGATGATGATG"),
    template = template
  )
  built <- plotly::plotly_build(p)

  expect_equal(unname(p$x$layout$xaxis$range), c(1, 12))
  expect_equal(vapply(p$x$layout$shapes, function(shape) shape$x1, numeric(1)), rep(12, 3))
  expect_equal(vapply(template$x$layout$shapes, function(shape) shape$x1, numeric(1)), rep(1, 3))
  expect_true(length(built$x$data) > length(template$x$data))
})

test_that("geneModelPanelPlotly exposes gene id labels as a named legend item", {
  dt <- data.table::data.table(
    gene_names = c("txA", "txA", "txB"),
    rect_starts = c(1, 8, 3),
    rect_ends = c(5, 12, 10),
    labels_locations = c(3, 10, 6.5),
    layers = c(1, 1, 2),
    type = c("cds", "cds", "cds"),
    cols = c("#F8766D", "#F8766D", "#00BA38")
  )

  built <- plotly::plotly_build(RiboCrypt:::geneModelPanelPlotly(dt))
  legend_traces <- Filter(function(tr) isTRUE(tr$showlegend), built$x$data)

  expect_length(legend_traces, 1)
  expect_identical(legend_traces[[1]]$name, "id")
  expect_identical(legend_traces[[1]]$legendgroup, "id")
  expect_identical(legend_traces[[1]]$mode, "text")
  expect_equal(sort(as.character(legend_traces[[1]]$text)), c("txA", "txB"))
})

test_that("browser_legend_cleanup keeps gene id legend item alongside frame legends", {
  profile <- data.table::data.table(
    position = 1:6,
    count = c(1, 2, 3, 2, 1, 0),
    frame = factor(c(0, 1, 2, 0, 1, 2))
  )
  gene_dt <- data.table::data.table(
    gene_names = c("txA", "txA"),
    rect_starts = c(1, 8),
    rect_ends = c(5, 12),
    labels_locations = c(3, 10),
    layers = c(1, 1),
    type = c("cds", "cds"),
    cols = c("#F8766D", "#F8766D")
  )

  p1 <- RiboCrypt:::createSinglePlot(
    profile, TRUE, "R", NULL, "a", "a", numeric(),
    type = "lines", lib_index = 1, total_libs = 1
  )
  gene_plot <- RiboCrypt:::geneModelPanelPlotly(gene_dt)

  polished <- suppressWarnings(RiboCrypt:::browser_legend_cleanup(
    plotly::subplot(list(p1, gene_plot), nrows = 2, shareX = TRUE, titleY = TRUE, titleX = TRUE)
  ))

  legend_traces <- Filter(function(tr) isTRUE(tr$showlegend), polished$x$data)
  legend_names <- vapply(legend_traces, function(tr) as.character(if (is.null(tr$name)) "" else tr$name), character(1))

  expect_equal(sort(legend_names[nzchar(legend_names)]), c("0", "1", "2", "id"))
})

test_that("ntSeqPanelPlotly reuses template x-axis layout without mutating template", {
  template <- RiboCrypt:::ntSeqPanelPlotlyTemplate()

  p <- RiboCrypt:::ntSeqPanelPlotly(Biostrings::DNAString("ATGATG"), template = template)

  expect_equal(unname(p$x$layout$xaxis$range), c(1, 6))
  expect_equal(unname(template$x$layout$xaxis$range), c(1, 1))
})

test_that("multiOmicsPlot_bottom_panels uses nt sequence plotly template", {
  display_range <- tx[1]
  template <- RiboCrypt:::ntSeqPanelPlotlyTemplate()
  aa_template <- RiboCrypt:::aaSeqPanelPlotlyTemplate()

  panels <- RiboCrypt:::multiOmicsPlot_bottom_panels(
    reference_sequence = ORFik::findFa(df),
    display_range = display_range,
    annotation = cds,
    start_codons = "ATG",
    stop_codons = c("TAA", "TAG", "TGA"),
    custom_motif = NULL,
    custom_regions = NULL,
    viewMode = "tx",
    tx_annotation = tx[1],
    collapse_intron_flank = 100,
    frame_colors = "R",
    templates = list(
      nt_seq_panel_plotly = template,
      aa_seq_panel_plotly = aa_template
    )
  )

  expect_equal(
    unname(panels$seq_nt_panel$x$layout$xaxis$range),
    c(1, Biostrings::nchar(panels$target_seq[[1]]))
  )
  expect_equal(unname(template$x$layout$xaxis$range), c(1, 1))
  expect_equal(
    vapply(panels$seq_panel$x$layout$shapes, function(shape) shape$x1, numeric(1)),
    rep(Biostrings::nchar(panels$target_seq[[1]]), 3)
  )
  expect_equal(vapply(aa_template$x$layout$shapes, function(shape) shape$x1, numeric(1)), rep(1, 3))
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
  expect_identical(polished$x$layout$legend$x, 1.02)
  expect_identical(polished$x$layout$legend$xanchor, "left")
  expect_identical(polished$x$layout$legend$y, 0.92)
  expect_identical(polished$x$layout$legend$yanchor, "top")
  expect_identical(polished$x$layout$legend$orientation, "v")
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

test_that("browser_legend_cleanup keeps one legend item per shared frame", {
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

  polished <- RiboCrypt:::browser_legend_cleanup(
    plotly::subplot(list(p1, p2), nrows = 2, shareX = TRUE, titleY = TRUE, titleX = TRUE)
  )

  legend_traces <- Filter(function(tr) isTRUE(tr$showlegend), polished$x$data)
  legend_names <- vapply(legend_traces, function(tr) as.character(if (is.null(tr$name)) "" else tr$name), character(1))

  expect_equal(sort(legend_names[nzchar(legend_names)]), c("0", "1", "2"))
})

test_that("createSinglePlot uses native plotly traces for supported track types", {
  profile <- data.table::data.table(
    position = rep(1:4, each = 3),
    count = c(1, 2, 3, 2, 1, 2, 3, 2, 1, 1, 3, 2),
    frame = factor(rep(0:2, 4))
  )

  line_plot <- plotly::plotly_build(RiboCrypt:::createSinglePlot(
    profile, TRUE, "R", NULL, "a", "a", numeric(),
    type = "lines", lib_index = 1, total_libs = 2
  ))
  line_traces <- Filter(function(tr) isTRUE(tr$showlegend), line_plot$x$data)
  expect_true(length(line_traces) >= 3)
  expect_true(all(vapply(line_traces, function(tr) identical(tr$type, "scatter"), logical(1))))
  expect_true(all(vapply(line_traces, function(tr) identical(tr$mode, "lines"), logical(1))))
  expect_true(all(vapply(line_traces, function(tr) any(grepl("%\\{x:\\.0f\\}", tr$hovertemplate)), logical(1))))

  column_plot <- plotly::plotly_build(RiboCrypt:::createSinglePlot(
    profile, TRUE, "R", NULL, "a", "a", numeric(),
    type = "columns", lib_index = 1, total_libs = 2
  ))
  column_traces <- Filter(function(tr) isTRUE(tr$showlegend), column_plot$x$data)
  expect_true(all(vapply(column_traces, function(tr) identical(tr$type, "scattergl"), logical(1))))
  expect_true(all(vapply(column_traces, function(tr) identical(tr$mode, "lines"), logical(1))))
  expect_true(all(vapply(column_traces, function(tr) identical(tr$line$simplify, FALSE), logical(1))))
  expect_true(all(vapply(column_traces, function(tr) identical(tr$visible, "legendonly"), logical(1))))
  expect_true(all(vapply(column_traces, function(tr) identical(tr$meta$rc_columns_switch, "gl_subset"), logical(1))))
  hidden_column_line_traces <- Filter(function(tr) identical(tr$meta$rc_columns_switch, "line"), column_plot$x$data)
  expect_length(hidden_column_line_traces, 3)
  expect_true(all(vapply(hidden_column_line_traces, function(tr) identical(tr$visible, TRUE), logical(1))))
  expect_equal(vapply(hidden_column_line_traces, function(tr) length(tr$x), integer(1)), c(4L, 4L, 4L))

  column_gl_plot <- plotly::plotly_build(RiboCrypt:::createSinglePlot(
    profile, TRUE, "R", NULL, "a", "a", numeric(),
    type = "columns", lib_index = 1, total_libs = 2,
    templates = list(cov_panel_columns_plotly = RiboCrypt:::covPanelColumnsGLPlotlyTemplate())
  ))
  column_gl_traces <- Filter(function(tr) isTRUE(tr$showlegend), column_gl_plot$x$data)
  expect_true(all(vapply(column_gl_traces, function(tr) identical(tr$type, "scattergl"), logical(1))))
  expect_true(all(vapply(column_gl_traces, function(tr) identical(tr$mode, "lines"), logical(1))))
  expect_true(all(vapply(column_gl_traces, function(tr) identical(tr$line$simplify, FALSE), logical(1))))

  area_plot <- plotly::plotly_build(RiboCrypt:::createSinglePlot(
    profile, TRUE, "R", NULL, "a", "a", numeric(),
    type = "area", lib_index = 1, total_libs = 2
  ))
  area_traces <- Filter(function(tr) isTRUE(tr$showlegend), area_plot$x$data)
  expect_true(all(vapply(area_traces, function(tr) identical(tr$type, "scatter"), logical(1))))
  expect_true(all(vapply(area_traces, function(tr) identical(tr$fill, "tozeroy"), logical(1))))
  expect_true(all(vapply(area_traces, function(tr) {
    if (!grepl("^#[0-9A-Fa-f]{8}$", tr$fillcolor)) return(FALSE)
    alpha_hex <- substr(tr$fillcolor, 8, 9)
    alpha <- strtoi(alpha_hex, base = 16) / 255
    alpha >= 0.58 && alpha <= 0.62
  }, logical(1))))

  stack_plot <- plotly::plotly_build(RiboCrypt:::createSinglePlot(
    profile, TRUE, "R", NULL, "a", "a", numeric(),
    type = "stacks", lib_index = 1, total_libs = 2
  ))
  stack_traces <- Filter(function(tr) isTRUE(tr$showlegend), stack_plot$x$data)
  expect_true(all(vapply(stack_traces, function(tr) identical(tr$type, "scatter"), logical(1))))
  expect_true(all(vapply(stack_traces, function(tr) identical(tr$stackgroup, "coverage"), logical(1))))

  heatmap_plot <- plotly::plotly_build(RiboCrypt:::createSinglePlot(
    profile, TRUE, "R", NULL, "a", "a", numeric(),
    type = "heatmap", lib_index = 1, total_libs = 2
  ))
  expect_identical(heatmap_plot$x$data[[1]]$type, "heatmap")
  expect_false(isTRUE(heatmap_plot$x$data[[1]]$showscale))
})

test_that("multiOmicsPlot_all_track_plots assigns distinct lib_index groups for multiple column tracks", {
  profiles <- list(
    data.table::data.table(
      position = rep(1:4, each = 3),
      count = c(1, 2, 3, 2, 1, 2, 3, 2, 1, 1, 3, 2),
      frame = factor(rep(0:2, 4))
    ),
    data.table::data.table(
      position = rep(1:4, each = 3),
      count = c(2, 1, 2, 3, 2, 1, 1, 3, 2, 2, 2, 1),
      frame = factor(rep(0:2, 4))
    )
  )

  track_panel <- RiboCrypt:::multiOmicsPlot_all_track_plots(
    profiles = profiles,
    withFrames = c(TRUE, TRUE),
    frame_colors = "R",
    colors = c(NA, NA),
    ylabels = c("lib1", "lib2"),
    ylabels_full_name = c("lib1", "lib2"),
    lines = numeric(),
    frames_type = "columns",
    summary_track = FALSE,
    summary_track_type = "columns",
    BPPARAM = BiocParallel::SerialParam(),
    templates = list(cov_panel_columns_plotly = RiboCrypt:::covPanelColumnsPlotlyTemplate())
  )

  expect_length(track_panel$plots, 2)
  group_ids <- vapply(track_panel$plots, function(plot) plot$x$data[[1]]$meta$rc_columns_group, character(1))
  expect_equal(group_ids, c("track_1", "track_2"))
})

test_that("createSinglePlot reuses framed coverage template for line tracks", {
  profile <- data.table::data.table(
    position = rep(1:4, each = 3),
    count = c(1, 2, 3, 2, 1, 2, 3, 2, 1, 1, 3, 2),
    frame = factor(rep(0:2, 4))
  )
  template <- RiboCrypt:::covPanelWithFramesPlotlyTemplate()

  p <- RiboCrypt:::createSinglePlot(
    profile, TRUE, "R", NULL, "a", "a", numeric(),
    type = "lines", lib_index = 1, total_libs = 2,
    templates = list(cov_panel_with_frames_plotly = template)
  )
  built <- plotly::plotly_build(p)

  expect_equal(length(built$x$data), 4)
  expect_equal(vapply(built$x$data[1:3], function(tr) length(tr$x), integer(1)), c(4L, 4L, 4L))
  expect_equal(vapply(template$x$data, function(tr) length(tr$x), integer(1)), c(1L, 1L, 1L))
})

test_that("createSinglePlot reuses non-framed coverage template for line tracks", {
  profile <- data.table::data.table(
    position = 1:4,
    count = c(1, 2, 3, 2)
  )
  template <- RiboCrypt:::covPanelWithoutFramesPlotlyTemplate()

  p <- RiboCrypt:::createSinglePlot(
    profile, FALSE, "R", "#4C78A8", "a", "a", numeric(),
    type = "lines", lib_index = 1, total_libs = 2,
    templates = list(cov_panel_without_frames_plotly = template)
  )
  built <- plotly::plotly_build(p)

  expect_true(length(built$x$data) >= 2)
  expect_equal(as.numeric(built$x$data[[1]]$x), 1:4)
  expect_equal(as.numeric(built$x$data[[1]]$y), c(1, 2, 3, 2))
  expect_identical(length(template$x$data[[1]]$x), 1L)
})

test_that("createSinglePlot reuses column coverage template", {
  profile <- data.table::data.table(
    position = rep(1:4, each = 3),
    count = c(1, 2, 3, 2, 1, 2, 3, 2, 1, 1, 3, 2),
    frame = factor(rep(0:2, 4))
  )
  template <- RiboCrypt:::covPanelColumnsPlotlyTemplate()

  p <- RiboCrypt:::createSinglePlot(
    profile, TRUE, "R", NULL, "a", "a", numeric(),
    type = "columns", lib_index = 1, total_libs = 2,
    templates = list(cov_panel_columns_plotly = template)
  )
  built <- plotly::plotly_build(p)

  expect_true(all(vapply(built$x$data[1:3], function(tr) identical(tr$type, "scattergl"), logical(1))))
  expect_true(all(vapply(built$x$data[1:3], function(tr) identical(tr$visible, "legendonly"), logical(1))))
  expect_equal(vapply(built$x$data[1:3], function(tr) length(tr$x), integer(1)), c(1L, 1L, 1L))
  expect_true(all(vapply(built$x$data[4:6], function(tr) identical(tr$type, "scatter"), logical(1))))
  expect_true(all(vapply(built$x$data[4:6], function(tr) identical(tr$visible, TRUE), logical(1))))
  expect_equal(vapply(built$x$data[4:6], function(tr) length(tr$x), integer(1)), c(4L, 4L, 4L))
  expect_equal(vapply(template$x$data, function(tr) length(tr$x), integer(1)), rep(1L, 6))
})

test_that("createSinglePlot reuses GL column coverage template", {
  profile <- data.table::data.table(
    position = rep(1:4, each = 3),
    count = c(1, 2, 3, 2, 1, 2, 3, 2, 1, 1, 3, 2),
    frame = factor(rep(0:2, 4))
  )
  template <- RiboCrypt:::covPanelColumnsGLPlotlyTemplate()

  p <- RiboCrypt:::createSinglePlot(
    profile, TRUE, "R", NULL, "a", "a", numeric(),
    type = "columns", lib_index = 1, total_libs = 2,
    templates = list(cov_panel_columns_plotly = template)
  )
  built <- plotly::plotly_build(p)

  expect_true(all(vapply(built$x$data[1:3], function(tr) identical(tr$type, "scattergl"), logical(1))))
  expect_true(all(vapply(built$x$data[1:3], function(tr) identical(tr$mode, "lines"), logical(1))))
  expect_equal(vapply(template$x$data, function(tr) length(tr$x), integer(1)), c(1L, 1L, 1L))
})

test_that("browser_plot_final_layout_polish adds a render hook for column zoom switching", {
  p <- plotly::subplot(
    list(
      RiboCrypt:::createSinglePlot(
        data.table::data.table(
          position = rep(1:4, each = 3),
          count = c(1, 2, 3, 2, 1, 2, 3, 2, 1, 1, 3, 2),
          frame = factor(rep(0:2, 4))
        ),
        TRUE, "R", NULL, "a", "a", numeric(),
        type = "columns", lib_index = 1, total_libs = 1
      )
    ),
    nrows = 1, shareX = TRUE, titleY = TRUE, titleX = TRUE
  )

  polished <- RiboCrypt:::browser_plot_final_layout_polish(
    p, "default", tx[1], NULL, NULL, "svg", NULL, numeric(0), 1
  )

  render_hooks <- polished$jsHooks$render
  expect_true(length(render_hooks) >= 2)
  render_code <- paste(vapply(render_hooks, `[[`, character(1), "code"), collapse = "\n")
  expect_match(render_code, "gl_subset")
  expect_match(render_code, "buildSegments")
  expect_match(render_code, "getGlLineWidth")
  expect_match(render_code, "line\\.width")
  expect_match(render_code, "hoverinfo")
  expect_match(render_code, "hovertemplate")
  expect_match(render_code, "%\\{text\\}<extra></extra>")
  expect_match(render_code, "max_x")
  expect_match(render_code, "clampRange")
})


test_that("createSinglePlot reuses area coverage template", {
  profile <- data.table::data.table(
    position = rep(1:4, each = 3),
    count = c(1, 2, 3, 2, 1, 2, 3, 2, 1, 1, 3, 2),
    frame = factor(rep(0:2, 4))
  )
  template <- RiboCrypt:::covPanelAreaPlotlyTemplate()

  p <- RiboCrypt:::createSinglePlot(
    profile, TRUE, "R", NULL, "a", "a", numeric(),
    type = "area", lib_index = 1, total_libs = 2,
    templates = list(cov_panel_area_plotly = template)
  )
  built <- plotly::plotly_build(p)

  expect_true(all(vapply(built$x$data[1:3], function(tr) identical(tr$fill, "tozeroy"), logical(1))))
  expect_equal(vapply(built$x$data[1:3], function(tr) length(tr$x), integer(1)), c(4L, 4L, 4L))
  expect_equal(vapply(template$x$data, function(tr) length(tr$x), integer(1)), c(1L, 1L, 1L))
})

test_that("createSinglePlot reuses heatmap coverage template", {
  profile <- data.table::data.table(
    position = rep(1:4, each = 3),
    count = c(1, 2, 3, 2, 1, 2, 3, 2, 1, 1, 3, 2),
    frame = factor(rep(0:2, 4))
  )
  template <- RiboCrypt:::covPanelHeatmapPlotlyTemplate()

  p <- RiboCrypt:::createSinglePlot(
    profile, TRUE, "R", NULL, "a", "a", numeric(),
    type = "heatmap", lib_index = 1, total_libs = 2,
    templates = list(cov_panel_heatmap_plotly = template)
  )
  built <- plotly::plotly_build(p)

  expect_identical(built$x$data[[1]]$type, "heatmap")
  expect_equal(as.numeric(built$x$data[[1]]$x), rep(1:4, each = 3))
  expect_identical(dim(built$x$data[[1]]$z), c(1L, 12L))
  expect_identical(length(template$x$data[[1]]$x), 1L)
})

test_that("getPlotAnimate produces native plotly animation frames", {
  profile <- data.table::data.table(
    position = rep(1:4, each = 3),
    count = c(1, 2, 3, 2, 1, 2, 3, 2, 1, 1, 3, 2),
    frame = factor(rep(0:2, 4))
  )
  profile_anim <- data.table::rbindlist(list(
    a = profile,
    b = data.table::copy(profile)[, count := count + 1L]
  ), idcol = "file")

  animate_plot <- suppressWarnings(plotly::plotly_build(
    RiboCrypt:::getPlotAnimate(profile_anim, TRUE, NULL, "R", "anim", numeric())
  ))

  expect_s3_class(animate_plot, "plotly")
  expect_equal(length(animate_plot$x$frames), 2)
  expect_true(all(vapply(animate_plot$x$frames[[1]]$data, function(tr) identical(tr$type, "scatter"), logical(1))))
  expect_true(all(vapply(animate_plot$x$frames[[1]]$data, function(tr) any(grepl("%\\{x:\\.0f\\}", tr$hovertemplate)), logical(1))))
})

test_that("lineDeSimplify does not warn on animated plotly tracks", {
  profile <- data.table::data.table(
    position = rep(1:4, each = 3),
    count = c(1, 2, 3, 2, 1, 2, 3, 2, 1, 1, 3, 2),
    frame = factor(rep(0:2, 4))
  )
  profile_anim <- data.table::rbindlist(list(
    a = profile,
    b = data.table::copy(profile)[, count := count + 1L]
  ), idcol = "file")

  expect_no_warning(
    RiboCrypt:::lineDeSimplify(
      RiboCrypt:::getPlotAnimate(profile_anim, TRUE, NULL, "R", "anim", numeric())
    )
  )
})
