test_that("Y-axis ranges accept auto, maxima and explicit limits", {
  expect_null(browser_parse_y_range(NULL))
  expect_null(browser_parse_y_range(" Auto "))
  for (value in c("500000", "500k", "0:500K", "0 : 5e5", "0:.5m")) {
    expect_identical(browser_parse_y_range(value), c(0, 500000))
  }
  expect_equal(browser_parse_y_range("1.5:2.5"), c(1.5, 2.5))
  for (value in c("", "bogus", "1:2:3", "1:", ":5", "NaN", "Inf", "0:Inf",
                  "-1:5", "5:-1", "0", "5:5", "10:5", "1e999", "500,000", "1;system('x')")) {
    expect_error(browser_parse_y_range(value), "Y-axis")
  }
  expect_error(browser_parse_y_range(NA), "single value")
  expect_error(browser_parse_y_range(c("1", "2")), "single value")
})

test_that("fixed ranges survive Plotly build without touching non-coverage axes", {
  make_plot <- function() plotly::plot_ly(x = 1:3, y = c(0, 2, 4), type = "scatter", mode = "lines")
  plot <- fast_subplot_shared_x(list(make_plot(), make_plot(), make_plot()), shareX = TRUE, nrows = 3)
  original <- unserialize(serialize(plot, NULL))
  fixed <- browser_apply_y_range(plot, c(0, 500000), c("area", "heatmap"))
  built <- suppressWarnings(plotly::plotly_build(fixed))
  expect_equal(built$x$layout$yaxis$range, c(0, 500000))
  expect_false(built$x$layout$yaxis$autorange)
  expect_equal(built$x$layout$yaxis$tickmode, "auto")
  expect_identical(fixed$x$layout$yaxis2, plot$x$layout$yaxis2)
  expect_identical(fixed$x$layout$yaxis3, plot$x$layout$yaxis3)
  expect_identical(browser_apply_y_range(plot, NULL, "area"), plot)
  expect_identical(plot, original)
  all_fixed <- browser_apply_y_range(plot, c(10, 20), c("lines", "columns", "stacks"))
  expect_true(all(vapply(all_fixed$x$layout[c("yaxis", "yaxis2", "yaxis3")],
                        function(axis) identical(axis$range, c(10, 20)), logical(1))))
})

test_that("summary and animation axes are counted in display order", {
  expect_equal(browser_coverage_track_types(list(frames_type = "area"), list(1, 2)), rep("area", 2))
  expect_equal(browser_coverage_track_types(list(frames_type = "animate", summary_track = TRUE,
                                                summary_track_type = "lines"), list(1, 2)), c("lines", "animate"))
  expect_equal(browser_coverage_track_types(list(frames_type = "heatmap", summary_track = TRUE,
                                                summary_track_type = "area"), list(1, 2)), c("area", "heatmap", "heatmap"))
})

test_that("Y-axis ranges invalidate coverage caches but not annotation caches", {
  df <- ORFik::ORFik.template.experiment()[9:10, ]
  auto <- hash_strings_browser(list(), df, ciw = 0)
  fixed <- hash_strings_browser(list(y_range = "500k"), df, ciw = 0)
  expect_identical(auto$hash_bottom, fixed$hash_bottom)
  expect_false(identical(auto$hash_browser, fixed$hash_browser))
  expect_identical(auto, hash_strings_browser(list(y_range = "auto"), df, ciw = 0))
})

test_that("both URL systems retain manual ranges and Observatory rejects invalid ones", {
  url <- make_url_from_inputs_parameters(list(y_range = "0:500k"))
  expect_match(url, "y_range=0%3A500k", fixed = TRUE)
  settings <- observatory_capture_browser_settings(list(y_range = "10:500k"))
  restored <- observatory_normalize_browser_settings(settings)
  expect_equal(restored$y_range, "10:500k")
  expect_true(observatory_browser_settings_ready(restored, list(y_range = "10:500k", frames_subset = character())))
  expect_false(observatory_browser_settings_ready(restored, list(y_range = "auto")))
  expect_error(observatory_normalize_browser_settings(list(y_range = "10:5")), "maximum greater")
})
