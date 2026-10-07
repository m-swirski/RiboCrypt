cell_details_fixture <- function() {
  raw <- matrix(c(0, 4, 8, 12, 2, 2, 0, 2), 4, 2, dimnames = list(NULL, c("a", "b")))
  table <- raw
  attr(table, "row_order_list") <- list(2, 1)
  attr(table, "km") <- list(cluster = c(1, 2))
  regions <- list(R1 = list(labels = "T1, TC1", kind = "Translon", ranges = IRanges::IRanges(1, 1)),
    R2 = list(labels = "CDS", kind = "CDS", ranges = IRanges::IRanges(2, 4)))
  view <- megabrowser_translon_display(table, raw, regions)
  meta <- data.table(Run = c("b", "a"), cluster = c(2, 1), grouping = c("B", "A"))
  list(view = view, meta = meta, display = list(table = view, meta = meta),
    metadata = data.table(Run = c("a", "b"), CONDITION = c("WT", "stress"), BioProject = c("P1", "P2")))
}

test_that("cell selection follows saved orders and rejects stale or invalid coordinates", {
  f <- cell_details_fixture()
  selected <- megabrowser_cell_selection(list(x = 1, y = 1), f$view, f$display, f$meta)
  expect_equal(selected, list(region = "R1", column = "b", runs = "b"))
  for (event in list(NULL, list(x = NA, y = 1), list(x = 999, y = 1), list(x = 1, y = -1)))
    expect_null(megabrowser_cell_selection(event, f$view, f$display, f$meta))
  focused <- megabrowser_focused_display(f$view, f$meta, "1")
  collapsed <- megabrowser_collapsed_display(focused$table, focused$meta)
  expect_equal(megabrowser_cell_selection(list(x = 2, y = 1), f$view, collapsed, f$meta)$runs, "a")
})

test_that("support uses summed base coverage and differs between library and cluster rows", {
  f <- cell_details_fixture()
  support <- megabrowser_region_support(f$view, f$display, f$meta, 5, 3)
  expect_equal(as.vector(support$count), c(0, 1, 0, 0))
  expect_equal(as.vector(support$low), c(TRUE, FALSE, TRUE, TRUE))
  collapsed <- megabrowser_collapsed_display(f$view, f$meta)
  expect_true(all(megabrowser_region_support(f$view, collapsed, f$meta, 5, 3)$low))
  expect_error(megabrowser_region_support(f$view, f$display, f$meta, -1, 3), "non-negative")
  expect_error(megabrowser_region_support(f$view, f$display, f$meta, 10, 1.5), "positive integer")
  before <- f$view
  p <- megabrowser_plotly_heatmap(f$view, c("white", "blue"))
  marked <- plotly::plotly_build(megabrowser_support_plot(p, f$view, f$display, f$meta,
    list(region_support = TRUE, region_min_coverage = 5, region_min_libraries = 3)))
  expect_identical(f$view, before)
  expect_equal(length(marked$x$data), 2)
  expect_equal(marked$x$data[[2]]$marker$symbol, "x")
  template <- plotly::plotly_build(p)
  templated <- megabrowser_plotly_heatmap(f$view, c("white", "blue"), template)
  marked_template <- plotly::plotly_build(megabrowser_support_plot(templated, f$view, f$display, f$meta,
    list(region_support = TRUE, region_min_coverage = 5, region_min_libraries = 3)))
  expect_equal(length(marked_template$x$data), 2)
  expect_equal(length(marked_template$x$data[[2]]$x), 3)
  expect_equal(marked_template$x$data[[1]]$z, plotly::plotly_build(templated)$x$data[[1]]$z)
})

test_that("cell ratios reuse raw densities and retain valid zero numerators", {
  f <- cell_details_fixture()
  selection <- list(region = "R1", column = "a", runs = "a")
  libraries <- megabrowser_cell_libraries(selection, f$view, f$meta, f$metadata)
  expect_equal(libraries$Run, c("a", "b"))
  expect_equal(libraries$CONDITION, c("WT", "stress"))
  ratios <- megabrowser_cell_ratios(libraries, f$view, "R2")
  expect_equal(ratios$Ratio, c(0, 1.5))
  expect_equal(ratios$Log2FC, c(-Inf, log2(1.5)))
  expect_equal(ratios$Reference_density, c(8, 4/3))
  expect_equal(ratios$Selected, c(TRUE, FALSE))
  expect_equal(ratios$Reference_coverage_sum, c(24, 4))
  expect_equal(ratios$Supported, c(FALSE, TRUE))
  supported <- megabrowser_cell_ratios(libraries, f$view, "R2", minimum = 3)
  expect_equal(supported$Supported, c(FALSE, FALSE))
  expect_equal(supported$Ratio, ratios$Ratio)
  expect_error(megabrowser_cell_ratios(libraries, f$view, "R2", minimum = -1), "non-negative")
  f$view["R2", "b"] <- 0
  expect_true(is.na(megabrowser_cell_ratios(libraries, f$view, "R2")$Ratio[2]))
})

test_that("ratio metadata retains missing categories, zero ratios and study concentration", {
  data <- data.table(Run = letters[1:6], Selected = c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE),
    CONDITION = c("WT", "WT", "stress", "stress", NA, ""), BioProject = c("P1", "P1", "P2", "P3", NA, NA),
    Ratio = c(0, 1, 3, 4, NA, 0))
  summary <- megabrowser_ratio_metadata(data, "CONDITION")
  expect_equal(summary[Term == "WT", Median_log2FC], -1)
  expect_equal(summary[Term == "(Missing)", Median_log2FC], -Inf)
  expect_equal(summary[Term == "stress", Studies], 2L)
  expect_equal(summary[Term == "WT", `Largest study (%)`], 100)
  expect_equal(summary[Term == "(Missing)", Undefined], 1L)
  expect_equal(summary[Term == "(Missing)", Zero], 1L)
  expect_equal(sum(summary$Libraries), 6L)
  expect_equal(summary$BH, p.adjust(summary$P, "BH"))
})

test_that("region/reference distributions use log2 scores without pseudocounts", {
  data <- data.table(Selected = TRUE, Density = c(0, .5, 1, 2, 8000, NA),
    Ratio = c(0, .5, 1, 2, 8000, NA))
  plot <- plotly::plotly_build(megabrowser_cell_distribution(data, "Ratio"))
  expect_equal(as.numeric(plot$x$data[[1]]$y), c(-1, 0, 1, log2(8000)))
  expect_match(plot$x$layout$yaxis$title, "log2", fixed = TRUE)
  raw <- plotly::plotly_build(megabrowser_cell_distribution(data, "Density"))
  expect_equal(as.numeric(raw$x$data[[1]]$y), c(0, .5, 1, 2, 8000))
  expect_error(megabrowser_cell_distribution(data[1], "Ratio"), "positive region")
})

test_that("biological labels preserve custom names and distinguish region positions", {
  regions <- list(list(kind = "Translon", ranges = IRanges::IRanges(1, 6)),
    list(kind = "Translon", ranges = IRanges::IRanges(4, 12)),
    list(kind = "CDS", ranges = IRanges::IRanges(10, 30)),
    list(kind = "Translon", ranges = IRanges::IRanges(15, 18)),
    list(kind = "Translon", ranges = IRanges::IRanges(35, 39)),
    megabrowser_user_region("40", "Custom signal", 40))
  expect_equal(megabrowser_biological_labels(regions), c("uORF1", "uORF2", "CDS", "internal ORF1", "dORF1", "Custom signal"))
})

test_that("cell modal uses the session theme and restores the surrounding context", {
  f <- cell_details_fixture()
  previous <- shiny::getShinyOption("bootstrapTheme")
  session <- list(ns = identity, getCurrentTheme = function() bslib::bs_theme(version = 5),
    setCurrentTheme = function(theme) invisible(theme),
    sendModal = function(type, message) {
      expect_equal(type, "show")
      expect_match(message$html, 'data-bs-toggle="tab"')
      expect_equal(bslib::theme_version(shiny::getShinyOption("bootstrapTheme")), "5")
    })
  megabrowser_show_cell_modal(session, list(region = "R1", column = "a"),
    attr(f$view, "translon_regions"), "CONDITION", "CONDITION")
  expect_identical(shiny::getShinyOption("bootstrapTheme"), previous)
  current <- NULL
  session$getCurrentTheme <- function() current
  session$setCurrentTheme <- function(theme) current <<- theme
  megabrowser_show_cell_modal(session, list(region = "R1", column = "a"),
    attr(f$view, "translon_regions"), "CONDITION", "CONDITION")
  expect_null(current)
  expect_identical(shiny::getShinyOption("bootstrapTheme"), previous)
})

test_that("cell modal is lazy, ratios respond to fields, and display changes clear selection", {
  f <- cell_details_fixture()
  shiny::testServer(function(input, output, session) {
    view <- reactive(f$view)
    display <- reactive({input$version; f$display})
    selected <- megabrowser_cell_outputs(input, output, session, view, display, reactive(list(meta = f$meta)), f$metadata)
  }, {
    session$setInputs(version = 1)
    expect_null(selected())
    session$setInputs(`plotly_click-mb_mid` = '{"x":1,"y":1}', region_min_coverage = 1, enrichment_metadata = "grouping")
    expect_equal(selected()$runs, "b")
    session$setInputs(cell_metric = "Ratio", cell_reference = "R2", cell_metadata = "CONDITION", cell_tabs = "Ratio by metadata")
    expect_match(output$cell_summary$html, "1 libraries")
    widget <- jsonlite::fromJSON(output$cell_ratio_metadata)
    expect_true(widget$x$options$serverSide)
    expect_match(widget$x$container, "Largest study")
    session$setInputs(version = 2)
    expect_null(selected())
  })
})
