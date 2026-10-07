translon_panel_fixture <- function() {
  data.table(gene_names = c("CDS", "CDS", "T1", "TC1", "T2", "T3"), type = "cds",
    rect_starts = c(5L, 12L, 2L, 2L, 1L, 6L), rect_ends = c(9L, 15L, 7L, 7L, 3L, 8L))
}

test_that("Heatmap gene annotation includes user regions on both strands across exons", {
  for (orientation in c("+", "-")) {
    tx <- GenomicRanges::GRangesList(tx = GenomicRanges::GRanges("chr1",
      IRanges::IRanges(c(100, 200), c(109, 219)), orientation))
    cds <- GenomicRanges::GRangesList(tx = GenomicRanges::GRanges("chr1",
      IRanges::IRanges(202, 215), orientation))
    predicted <- GenomicRanges::GRangesList(T = GenomicRanges::GRanges("chr1",
      IRanges::IRanges(100, 105), orientation))
    controller <- list(display_region = tx, customRegions = predicted)
    custom <- list(megabrowser_user_region("8:23", "UR1", 30),
                   megabrowser_user_region("25", "Edited label", 30))
    merged <- megabrowser_annotation_regions(controller, custom)
    expect_identical(unname(merged$T), unname(predicted$T))
    expect_identical(names(merged), c("T", "UR1", "Edited label"))
    panel <- createGeneModelPanel(tx, cds, tx_annotation = tx, custom_regions = merged,
                                 viewMode = "tx", collapse_intron_flank = 0)[[1]]
    expect_equal(min(panel[gene_names == "UR1", rect_starts]), 8L)
    expect_equal(max(panel[gene_names == "UR1", rect_ends]), 23L)
    expect_equal(sum(panel[gene_names == "UR1", rect_ends - rect_starts + 1L]), 16L)
    expect_equal(panel[gene_names == "Edited label", rect_starts], 25L)
    plot <- plotly::plotly_build(geneModelPanelPlotly(panel))
    expect_match(paste(unlist(lapply(plot$x$data, `[[`, "text")), collapse = " "), "UR1")
    expect_identical(megabrowser_annotation_regions(controller, list()), predicted)
  }
})

test_that("user regions accept inclusive coordinates and editable labels", {
  region <- megabrowser_user_region("40:80", "My ORF", 100)
  expect_equal(region$ranges, IRanges::IRanges(40, 80))
  expect_equal(region$labels, "My ORF")
  expect_identical(megabrowser_user_region("40", "UR1", 100), megabrowser_user_region("40:40", "UR1", 100))
  expect_equal(megabrowser_user_region(" 40 : 80 ", " UR1 ", 100)$labels, "UR1")
  expect_error(megabrowser_user_region("40:39", "UR1", 100), "Start must not exceed")
  for (invalid in c("", "abc", "1.5", "1:2:3", "-1", "1e2"))
    expect_error(megabrowser_user_region(invalid, "UR1", 100), "Enter a position")
  for (invalid in c("0", "101", "1:101", "999999999999999999999999"))
    expect_error(megabrowser_user_region(invalid, "UR1", 100), "between 1 and 100")
  expect_error(megabrowser_user_region("40", "  ", 100), "region label")
})

test_that("custom regions join sequential IDs, raw ratios and exports", {
  panel <- data.table(gene_names = c("T", "CDS"), type = "cds", rect_starts = c(1L, 5L), rect_ends = c(2L, 10L))
  testthat::local_mocked_bindings(megabrowser_translon_panel = function(...) panel, .package = "RiboCrypt")
  raw <- matrix(2, 10, 2, dimnames = list(NULL, c("a", "b")))
  raw[3:4, ] <- 8
  result <- megabrowser_translon_analysis(list(annotation = list(CDS = NULL)), list(table = raw),
    list(meta = data.table(Run = c("a", "b"), cluster = c(1, 2), grouping = "x")),
    custom = list(megabrowser_user_region("3:4", "My region", 10)), raw = raw)
  expect_equal(result$regions$Region, c("R1", "R2", "R3"))
  expect_equal(result$regions[Region == "R3", Labels], "My region")
  expect_equal(nrow(result$pairs), 3L)
  expect_equal(result$ratios[Numerator_label == "My region", Ratio], c(4, 4))
  expect_equal(nrow(result$statistics), 6L)
})

test_that("region popup rejects errors, accepts edited labels and increments defaults", {
  shiny::testServer(function(input, output, session) {
    custom <- reactiveVal(list())
    megabrowser_user_region_events(input, session, custom, reactive(matrix(1, 100, 2)))
  }, {
    session$setInputs(translon_add_region = 1, translon_region_coordinates = "40:39", translon_region_label = "UR1")
    session$setInputs(translon_region_save = 1)
    expect_length(custom(), 0)
    expect_match(output$translon_region_error$html, "Start must not exceed")
    session$setInputs(translon_region_coordinates = "40:80", translon_region_label = "My ORF", translon_region_save = 2)
    expect_equal(custom()[[1]]$labels, "My ORF")
    session$setInputs(translon_add_region = 2, translon_region_coordinates = "40", translon_region_label = "UR2", translon_region_save = 3)
    expect_length(custom(), 2)
    expect_equal(custom()[[2]]$ranges, IRanges::IRanges(40, 40))
  })
})

test_that("annotation duplicates retain aliases and clean CDS removes overlapping uORF bases", {
  regions <- megabrowser_translon_regions(translon_panel_fixture(), "CDS")
  table <- megabrowser_translon_region_table(regions)
  expect_equal(sum(grepl("T1", table$Labels)), 1L)
  expect_match(table[grepl("T1", Labels), Labels], "TC1", fixed = TRUE)
  expect_equal(table[Type == "clean_cds", Intervals], "8-9;12-15")
  expect_equal(table[Type == "clean_cds", Length_nt], 6L)
  expect_equal(table[Type == "CDS", Length_nt], 9L)
  duplicated_cds <- rbind(translon_panel_fixture(), data.table(gene_names = "TC2", type = "cds", rect_starts = 5L, rect_ends = 9L))
  expect_equal(nrow(megabrowser_translon_region_table(megabrowser_translon_regions(duplicated_cds, "CDS"))), nrow(table) + 1L)
})

test_that("clean CDS is omitted without overlap and can be empty when entirely covered", {
  panel <- data.table(gene_names = c("CDS", "T"), type = "cds", rect_starts = c(5L, 1L), rect_ends = c(9L, 3L))
  expect_false(any(megabrowser_translon_region_table(megabrowser_translon_regions(panel, "CDS"))$Type == "clean_cds"))
  panel[gene_names == "T", rect_ends := 12L]
  regions <- megabrowser_translon_regions(panel, "CDS")
  expect_equal(megabrowser_translon_region_table(regions)[Type == "clean_cds", Length_nt], 0L)
  raw <- matrix(1, 15, 2, dimnames = list(NULL, c("a", "b")))
  expect_true(all(is.na(megabrowser_translon_density(raw, regions)[, names(regions)[vapply(regions, function(region) region$kind == "clean_cds", logical(1))]])))
})

test_that("minus-strand annotation track coordinates give the correct upstream overlap", {
  tx <- GenomicRanges::GRangesList(tx = GenomicRanges::GRanges("chr1", IRanges::IRanges(100, 199), "-"))
  cds <- GenomicRanges::GRangesList(tx = GenomicRanges::GRanges("chr1", IRanges::IRanges(120, 159), "-"))
  uorf <- GenomicRanges::GRanges("chr1", IRanges::IRanges(150, 179), "-")
  panel <- createGeneModelPanel(tx, cds, tx_annotation = tx,
    custom_regions = GenomicRanges::GRangesList(T = uorf, TC = uorf), viewMode = "tx", collapse_intron_flank = 0)[[1]]
  regions <- megabrowser_translon_region_table(megabrowser_translon_regions(panel, "tx"))
  expect_equal(regions[Type == "CDS", Intervals], "41-80")
  expect_equal(regions[Type == "clean_cds", Intervals], "51-80")
  expect_equal(regions[Type == "Translon", Length_nt], 30L)
  expect_equal(regions[Type == "Translon", Labels], "T, TC")
})

test_that("raw density uses exact exons and all pairs exclude zero denominators", {
  regions <- list(R1 = list(labels = "uORF", kind = "Translon", ranges = IRanges::IRanges(1, 2)),
    R2 = list(labels = "CDS", kind = "CDS", ranges = IRanges::IRanges(c(5, 9), c(6, 10))))
  raw <- matrix(100, 10, 3, dimnames = list(NULL, c("a", "b", "c")))
  raw[1:2, ] <- rep(c(4, 0, 2), each = 2)
  raw[c(5:6, 9:10), ] <- rep(c(2, 3, 0), each = 4)
  density <- megabrowser_translon_density(raw, regions)
  expect_equal(density[, "R2"], setNames(c(2, 3, 0), c("a", "b", "c")))
  pairs <- megabrowser_translon_pairs(regions)
  ratios <- megabrowser_translon_ratios(density, pairs, list(A = c(1, 3), B = 2))
  expect_equal(ratios$Ratio, c(2, 0, NA))
  expect_equal(ratios$Group, c("A", "B", "A"))
  expect_equal(ratios$Log2_ratio[1], 1)
  expect_true(all(is.na(ratios$Log2_ratio[2:3])))
  stats <- megabrowser_translon_stats(ratios)
  expect_equal(stats[Group == "A", Undefined], 1L)
  expect_equal(stats[Group == "B", Zero], 1L)
  expect_true(all(is.na(stats$P)))
})

test_that("rank statistics and BH adjustment compare full library groups against rest", {
  ratios <- data.table(Pair = "R1/R2", Group = rep(c("A", "B"), each = 3), Ratio = c(4, 5, 6, 1, 2, 3))
  stats <- megabrowser_translon_stats(ratios)
  expect_equal(stats$Median, c(5, 2))
  expect_equal(stats$Rank_biserial, c(1, -1))
  expect_equal(stats$P, rep(suppressWarnings(wilcox.test(4:6, 1:3, exact = FALSE)$p.value), 2))
  expect_equal(stats$BH, p.adjust(stats$P, "BH"))
})

test_that("analysis reloads raw coverage rather than taking ratios of heatmap logs", {
  raw <- matrix(0, 10, 2, dimnames = list(NULL, c("a", "b")))
  raw[1:2, ] <- rep(c(4, 8), each = 2)
  raw[5:10, ] <- rep(c(2, 2), each = 6)
  normalized <- matrix(99, 10, 2, dimnames = dimnames(raw))
  meta <- data.table(Run = c("b", "a"), cluster = c(2, 1), grouping = "x")
  panel <- data.table(gene_names = c("T", "CDS"), type = "cds", rect_starts = c(1L, 5L), rect_ends = c(2L, 10L))
  testthat::local_mocked_bindings(megabrowser_translon_panel = function(...) panel,
    load_collection = function(path, columns) {
      expect_equal(columns, c("a", "b"))
      raw
    }, .package = "RiboCrypt")
  result <- megabrowser_translon_analysis(list(annotation = setNames(list(NULL), "CDS"), table_path = "mock"),
    list(table = normalized), list(meta = meta))
  expect_equal(result$ratios$Ratio, c(2, 4))
  expect_equal(result$ratios$Group, c("1", "2"))
  expect_equal(normalized, matrix(99, 10, 2, dimnames = dimnames(raw)))
})

test_that("insufficient annotation gives useful validation rather than a blank page", {
  expect_error(megabrowser_translon_regions(data.table(type = "utr"), "CDS"), "No coding regions")
  expect_error(megabrowser_translon_pairs(list(R1 = list(kind = "CDS"))), "At least two distinct")
  expect_error(megabrowser_translon_ratio_plot(data.table(Log2_ratio = NA_real_)), "No positive, finite ratios")
})

test_that("translon analysis is gated by the tab and ignores collapse and group focus", {
  calls <- 0L
  result <- list(regions = data.table(Region = c("R1", "R2"), Labels = c("T", "CDS")),
    pairs = data.table(Pair = "R1/R2", Numerator = "R1", Denominator = "R2"), ratios = data.table(), statistics = data.table())
  testthat::local_mocked_bindings(megabrowser_translon_analysis = function(...) {calls <<- calls + 1L; result}, .package = "RiboCrypt")
  shiny::testServer(function(input, output, session) {
    analysis <- megabrowser_translon_outputs(input, output, session,
      reactive(list(table_hash = if (is.null(input$hash)) "translon-test" else input$hash, table_plot_hash = "plot")), reactive(NULL), reactive(NULL))
  }, {
    session$setInputs(mb_tabs = "Heatmap")
    expect_equal(calls, 0L)
    session$setInputs(mb_tabs = "Translon enrichment")
    expect_identical(analysis(), result)
    expect_equal(calls, 1L)
    session$setInputs(collapsed_clusters = TRUE, visible_groups = "A")
    expect_identical(analysis(), result)
    expect_equal(calls, 1L)
    session$setInputs(mb_tabs = "Heatmap")
    session$setInputs(hash = "new-transcript")
    expect_equal(calls, 1L)
    session$setInputs(mb_tabs = "Translon enrichment")
    expect_equal(calls, 2L)
  })
})
