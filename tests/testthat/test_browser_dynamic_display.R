test_that("dynamic thresholds are precise", {
  f <- RiboCrypt:::browser_dynamic_style
  expect_equal(f(1), list(type = "columns", kmer = 1L))
  expect_equal(f(500), list(type = "columns", kmer = 1L))
  expect_equal(f(501), list(type = "area", kmer = 9L))
  expect_equal(f(10000), list(type = "area", kmer = 9L))
  expect_equal(f(10001), list(type = "lines", kmer = 9L))
})

test_that("dynamic view retains zoom and rejects events from another region", {
  region <- GenomicRanges::GRangesList(tx = GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 12000), "+"))
  input <- list(display_mode = "dynamic", tx = "tx", extendLeaders = 0, extendTrailers = 0, frames_type = "heatmap", kmer = 7)
  settings <- RiboCrypt:::browser_dynamic_settings(input, region, numeric())
  expect_equal(settings$type, "lines")
  input$dynamic_view <- list(context = settings$context, range = c(100, 599))
  settings <- RiboCrypt:::browser_dynamic_settings(input, region, numeric())
  expect_equal(settings$type, "columns")
  expect_equal(settings$range, c(100, 599))
  input$tx <- "other"
  expect_equal(RiboCrypt:::browser_dynamic_settings(input, region, numeric())$type, "lines")
  input$display_mode <- "manual"
  expect_equal(RiboCrypt:::browser_dynamic_settings(input, region, numeric())$type, "heatmap")
  expect_equal(RiboCrypt:::browser_dynamic_settings(input, region, numeric())$kmer, 7)
})

test_that("old and new observatory URL settings retain their display intent", {
  f <- RiboCrypt:::observatory_normalize_browser_settings
  expect_identical(f(list(frames_type = "area", kmer = 9))$display_mode, "manual")
  expect_identical(f(list(display_mode = "dynamic", kmer = 9))$display_mode, "dynamic")
  expect_error(f(list(display_mode = "bad")))
})

test_that("a new plot request invalidates a previous zoom event", {
  region <- GenomicRanges::GRangesList(tx = GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 12000), "+"))
  input <- list(go = 1L, display_mode = "dynamic", tx = "tx", extendLeaders = 0, extendTrailers = 0)
  first <- RiboCrypt:::browser_dynamic_settings(input, region, numeric())
  input$dynamic_view <- list(context = first$context, range = c(1, 300))
  expect_equal(RiboCrypt:::browser_dynamic_settings(input, region, numeric())$type, "columns")
  input$go <- 2L
  expect_equal(RiboCrypt:::browser_dynamic_settings(input, region, numeric())$type, "lines")
})
