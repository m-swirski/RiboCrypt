test_that("automatic plots preload the widget's own unmodified dependencies", {
  expected <- plotly::plot_ly()$dependencies
  expect_identical(plotly_startup_dependencies(TRUE), expected)
  expect_identical(plotly_startup_dependencies("TRUE"), expected)
  expect_identical(plotly_startup_dependencies("true"), expected)
  expect_true("plotly-main" %in% vapply(expected, `[[`, "", "name"))
  for (disabled in list(FALSE, "FALSE", "false", NULL, NA)) {
    expect_null(plotly_startup_dependencies(disabled))
  }
})

test_that("startup dependencies are registered once with Shiny and preserve its jquery", {
  dependencies <- plotly_startup_dependencies(TRUE)
  ui <- htmltools::tagList(dependencies, dependencies,
    shiny::selectizeInput("gene", "Gene", "ATF4"),
    plotly::plotlyOutput("plot"), DT::DTOutput("table"))
  app <- cache_static_app_ui(shiny::shinyApp(ui, function(input, output, session) {}))
  req <- list(REQUEST_METHOD = "GET", PATH_INFO = "/", QUERY_STRING = "")
  response <- app$httpHandler(req)
  page <- xml2::read_html(response$content)
  scripts <- xml2::xml_attr(xml2::xml_find_all(page, ".//script[@src]"), "src")
  expect_equal(sum(grepl("plotly-main-", scripts, fixed = TRUE)), 1L)
  expect_equal(sum(grepl("/jquery.min.js", scripts, fixed = TRUE)), 1L)
  registry <- xml2::xml_text(xml2::xml_find_first(page,
    './/script[@type="application/html-dependencies"]'))
  main <- Filter(function(dep) dep$name == "plotly-main", dependencies)[[1]]
  expect_match(registry, paste0(main$name, "[", main$version, "]"), fixed = TRUE)
  expect_identical(app$httpHandler(req), response)
  req$QUERY_STRING <- "gene=ATF4&go=true"
  query_page <- xml2::read_html(app$httpHandler(req)$content)
  expect_identical(xml2::xml_attr(xml2::xml_find_all(query_page, ".//script[@src]"), "src"), scripts)
  lazy_ui <- htmltools::tagList(plotly_startup_dependencies(FALSE),
    shiny::selectizeInput("gene", "Gene", "ATF4"), plotly::plotlyOutput("plot"))
  lazy_app <- shiny::shinyApp(lazy_ui, function(input, output, session) {})
  lazy_page <- xml2::read_html(lazy_app$httpHandler(req)$content)
  lazy_scripts <- xml2::xml_attr(xml2::xml_find_all(lazy_page, ".//script[@src]"), "src")
  expect_false(any(grepl("plotly-main-", lazy_scripts, fixed = TRUE)))
  expect_identical(scripts[grepl("/jquery.min.js", scripts, fixed = TRUE)],
                   lazy_scripts[grepl("/jquery.min.js", lazy_scripts, fixed = TRUE)])
})
