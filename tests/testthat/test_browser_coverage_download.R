coverage_test_plot <- function(profiles, labels, summary = FALSE) {
  structure(list(), coverage_export = list(profiles = profiles, labels = labels, summary = summary))
}

test_that("exports retain transformed values, track names and summary ordering", {
  profiles <- list(data.table::data.table(position = 1:3, count = c(0, log2(3), 2)),
                   data.table::data.table(position = 1:3, count = c(0.5, 0, 1)))
  plot <- coverage_test_plot(profiles, c("Library A", "Group,B"))
  result <- browser_coverage_table(plot)
  expect_identical(names(result), c("position", "Library A", "Group,B"))
  expect_equal(result[[2]], profiles[[1]]$count)
  expect_equal(result[[3]], profiles[[2]]$count)
  summary <- browser_coverage_table(coverage_test_plot(profiles, c("A", "B"), TRUE))
  expect_identical(names(summary), c("position", "summary", "B", "A"))
  expect_equal(summary$summary, profiles[[1]]$count + profiles[[2]]$count)
  file <- tempfile(fileext = ".csv")
  on.exit(unlink(file))
  data.table::fwrite(result, file)
  expect_equal(data.table::fread(file), result)
})

test_that("the shared plot builder retains exactly its plotted profiles for both browsers", {
  captured <- NULL
  profiles <- list(data.table::data.table(position = 1:3, count = c(0, 1.5, 0)))
  testthat::local_mocked_bindings(
    multiOmicsPlot_all_profiles = function(...) profiles,
    multiOmicsPlot_all_track_plots = function(profiles, ...) {
      captured <<- profiles
      list()
    },
    multiOmicsPlot_complete_plot = function(...) plotly::plot_ly(),
    .package = "RiboCrypt"
  )
  for (observatory in c(FALSE, TRUE)) {
    controls <- list(reads = list("reads"), withFrames = TRUE, viewMode = FALSE,
                     frames_type = "lines", frame_colors = "R", kmerLength = 3,
                     is_cellphone = FALSE, log_scale = TRUE, summary_track = FALSE,
                     frames_subset = "red")
    result <- browser_track_panel_shiny(function() controls,
      list(display_range = NULL, annotation_layers = 1, ncustom = 0),
      session = NULL, ylabels = if (observatory) "Selected group" else "Library A",
      profiles = if (observatory) profiles else NULL)
    expect_identical(attr(result, "coverage_export")$profiles, captured)
    expect_equal(browser_coverage_table(result)[[2]], captured[[1]]$count)
    expect_equal(names(browser_coverage_table(result))[2],
                 if (observatory) "Selected group" else "Library A")
    expect_null(result$x$coverage_export)
  }
})

test_that("coverage downloads use the generated snapshot without startup extraction", {
  calls <- 0L
  controls <- list(display_region = list(`tx/1` = NULL))
  shiny::testServer(function(input, output, session) {
    plot <- shiny::eventReactive(input$go, {
      calls <<- calls + 1L
      coverage_test_plot(list(data.table::data.table(position = 1:2, count = c(0, input$count))), "group")
    })
    browser_coverage_download(output, function() controls, plot)
  }, {
    expect_equal(calls, 0L)
    session$setInputs(go = 1, count = 7)
    expect_equal(data.table::fread(output$download_coverage),
                 data.table::data.table(position = 1:2, group = c(0, 7)))
    expect_equal(calls, 1L)
    session$setInputs(count = 99)
    expect_equal(data.table::fread(output$download_coverage)$group, c(0, 7))
    session$setInputs(go = 2)
    expect_equal(data.table::fread(output$download_coverage)$group, c(0, 99))
  })
  expect_equal(browser_coverage_filename(controls), "RiboCrypt_tx_1_coverage.csv")
  expect_equal(browser_coverage_filename(list()), "RiboCrypt_region_coverage.csv")
})
