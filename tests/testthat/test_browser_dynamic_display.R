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

test_that("client sidecar retains exact frame-specific smoothed profiles", {
  raw <- data.table::data.table(position = 1:30, frame = factor(rep(0:2, 10)), library = factor(rep(1, 30)), count = as.numeric(1:30))
  smooth <- RiboCrypt:::smoothenMultiSampCoverage(data.table::copy(raw), 9L, split_by_frame = TRUE)
  config <- RiboCrypt:::browser_dynamic_client_config(list(raw), list(smooth), "R", "mode", "status")
  for (phase in 0:2) expect_equal(config$tracks[[1]][[as.character(phase)]]$y,
                                  smooth[as.character(frame) == as.character(phase)]$count)
  bad <- data.table::copy(smooth); bad$position[1] <- 99
  expect_error(RiboCrypt:::browser_dynamic_client_config(list(raw), list(bad), "R", "mode", "status"), "position mismatch")
})

test_that("downloads select the profile currently displayed by dynamic view", {
  raw <- data.table::data.table(position = 1:3, count = c(1,2,3))
  smooth <- data.table::data.table(position = 1:3, count = c(2,2,2))
  p <- list(); attr(p, "coverage_export") <- list(profiles = list(raw), smooth_profiles = list(smooth), labels = "sample", summary = FALSE)
  expect_equal(RiboCrypt:::browser_coverage_table(p, "columns")$sample, c(1,2,3))
  expect_equal(RiboCrypt:::browser_coverage_table(p, "area")$sample, c(2,2,2))
  expect_equal(RiboCrypt:::browser_coverage_table(p)$sample, c(1,2,3))
})

test_that("dynamic client uses the safe server path for unsupported track combinations", {
  f <- RiboCrypt:::browser_dynamic_client_supported
  expect_true(f(TRUE, c(TRUE, TRUE), FALSE))
  expect_false(f(TRUE, c(TRUE, FALSE), FALSE))
  expect_false(f(TRUE, TRUE, TRUE))
  expect_false(f(TRUE, logical(), FALSE))
  expect_false(f(TRUE, NA, FALSE))
  expect_false(f(FALSE, TRUE, FALSE))
})
