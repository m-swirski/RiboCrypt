test_that("wide startup prepares the same sequence placeholder before rendering", {
  plot <- plotly::plotly_build(plotly::plot_ly(x = 1:3, y = 1:3, type = "scatter", mode = "lines"))
  data <- list(sequence = paste(rep("A", 1000), collapse = ""),
               traces = list(list(distance = 250, yaxis = "y2")))
  prepared <- browser_initial_sequence_placeholder(plot, data)
  trace <- tail(prepared$x$data, 1)[[1]]
  expect_equal(as.numeric(trace$x), 500.5)
  expect_equal(as.numeric(trace$y), 0.5)
  expect_identical(trace$name, "sequence_placeholder")
  expect_identical(trace$yaxis, "y2")
  expect_identical(trace$type, "scatter")
  expect_identical(head(prepared$x$data, -1), plot$x$data)
  expect_identical(browser_initial_sequence_placeholder(plot, data, c(40, 80)), plot)
  expect_identical(browser_initial_sequence_placeholder(plot, data, c(1, 251)), plot)
  zoomed <- browser_initial_sequence_placeholder(plot, data, c(10.2, 500.3))
  expect_equal(as.numeric(tail(zoomed$x$data, 1)[[1]]$x), 255.5)
})

test_that("browser x-axis styles already match the sequence callback", {
  axis <- browser_style_xaxis(list(), "xaxis", "xaxis")
  expect_identical(axis$ticks, "outside")
  expect_true(axis$showticklabels)
  for (property in c("showline", "showgrid", "zeroline")) expect_false(axis[[property]])
  expect_identical(browser_style_xaxis(list(), "xaxis2", "xaxis")$ticks, "")
})
