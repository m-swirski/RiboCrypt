test_that("long SVG coverage traces simplify drawing without changing points", {
  make_profile <- function(n) data.table::data.table(
    position = seq_len(n), count = as.numeric(seq_len(n) %% 101),
    frame = factor((seq_len(n) - 1L) %% 3L))
  templates <- list(cov_panel_with_frames_plotly = RiboCrypt:::covPanelWithFramesPlotlyTemplate(),
                    cov_panel_area_plotly = RiboCrypt:::covPanelAreaPlotlyTemplate())
  for (mode in c("lines", "area")) for (n in c(15000L, 15003L)) {
    profile <- make_profile(n)
    plot <- RiboCrypt:::singlePlot_select_plot_type(profile, TRUE, "R", NULL,
      numeric(), mode, 1L, templates = templates)
    for (i in seq_len(3L)) {
      trace <- plot$x$data[[i]]
      expected <- profile[as.character(frame) == as.character(i-1L)]
      expect_identical(trace$x, expected$position)
      expect_identical(trace$y, expected$count)
      expect_identical(trace$line$simplify, nrow(expected) > 5000L)
    }
  }
})

test_that("short tracks and detailed column traces remain unsimplified", {
  expect_false(RiboCrypt:::browser_simplify_coverage(5000L))
  expect_true(RiboCrypt:::browser_simplify_coverage(5001L))
  profile <- data.table::data.table(position = 1:18000, count = 1,
    frame = factor((0:17999) %% 3L))
  plot <- RiboCrypt:::singlePlot_select_plot_type(profile, TRUE, "R", NULL,
    numeric(), "columns", 1L, templates = list(cov_panel_columns_plotly = RiboCrypt:::covPanelColumnsPlotlyTemplate()))
  expect_true(all(vapply(plot$x$data[1:6], function(t) identical(t$line$simplify, FALSE), logical(1))))
})


test_that("SVG simplification can be disabled for exact path rendering", {
  old <- options(RiboCrypt.browser_simplify = FALSE)
  on.exit(options(old), add = TRUE)
  expect_false(RiboCrypt:::browser_simplify_coverage(100000L))
  profile <- data.table::data.table(position = 1:3, count = 1:3, frame = factor(0:2))
  config <- RiboCrypt:::browser_dynamic_client_config(list(profile), list(profile), "R", "mode", "status")
  expect_identical(config$simplify, FALSE)
})

test_that("non-template and unframed paths retain complete long coverage", {
  n <- 15003L
  profile <- data.table::data.table(position = seq_len(n), count = as.numeric(seq_len(n) %% 37L),
    frame = factor((seq_len(n)-1L) %% 3L))
  for (framed in c(TRUE, FALSE)) for (mode in c("lines", "area")) {
    plot <- RiboCrypt:::singlePlot_select_plot_type(profile, framed, "R", NULL,
      numeric(), mode, 1L, templates = NULL)
    traces <- plotly::plotly_build(plot)$x$data
    count <- if (framed) 3L else 1L
    for (i in seq_len(count)) {
      expect_true(traces[[i]]$line$simplify)
      expected <- if (framed) profile[as.character(frame) == as.character(i-1L)] else profile
      expect_equal(as.numeric(traces[[i]]$x), as.numeric(expected$position))
      expect_equal(as.numeric(traces[[i]]$y), as.numeric(expected$count))
    }
  }
})
