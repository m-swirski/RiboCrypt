test_that("static UI caching preserves Shiny's complete response", {
  renders <- 0L
  ui <- htmltools::tagList(
    htmltools::singleton(htmltools::tags$head(htmltools::tags$script("window.once = true;"))),
    shiny::selectizeInput("gene", "Gene", "ATF4"),
    plotly::plotlyOutput("plot"), DT::DTOutput("table"))
  app <- shiny::shinyApp(ui, function(input, output, session) {})
  handler <- app$httpHandler
  app$httpHandler <- function(req) { renders <<- renders + 1L; handler(req) }
  app <- cache_static_app_ui(app)
  req <- list(REQUEST_METHOD = "GET", PATH_INFO = "/", QUERY_STRING = "")
  expected <- handler(req)
  first <- app$httpHandler(req)
  second <- app$httpHandler(req)
  expect_identical(first, expected)
  expect_identical(second, expected)
  expect_equal(renders, 1L)
  expect_match(first$content, "application/shiny-singletons", fixed = TRUE)
  expect_match(first$content, "application/html-dependencies", fixed = TRUE)
  expect_match(first$content, "window.once = true;", fixed = TRUE)
})

test_that("query strings and non-root requests retain the original handler", {
  calls <- 0L
  app <- cache_static_app_ui(list(httpHandler = function(req) {
    calls <<- calls + 1L
    list(status = 200L, content = req)
  }))
  plain <- list(REQUEST_METHOD = "GET", PATH_INFO = "/", QUERY_STRING = "")
  app$httpHandler(plain)
  queries <- c("gene=ATF4&go=true", "obs_state=abc", "_state_id_=abc", "showcase=1")
  for (query in queries) {
    req <- plain; req$QUERY_STRING <- query
    expect_identical(app$httpHandler(req)$content, req)
    expect_identical(app$httpHandler(req)$content, req)
  }
  for (path in c("/asset.js", "/session/test/dataobj")) {
    req <- plain; req$PATH_INFO <- path
    expect_identical(app$httpHandler(req)$content, req)
  }
  req <- plain; req$REQUEST_METHOD <- "POST"
  expect_identical(app$httpHandler(req)$content, req)
  expect_equal(calls, 12L)
  app$httpHandler(plain)
  expect_equal(calls, 12L)
})

test_that("UI responses are isolated by app instance and failures are not cached", {
  calls <- 0L
  original <- list(httpHandler = function(req) {
    calls <<- calls + 1L
    if (calls == 1L) return(NULL)
    list(status = if (calls == 2L) 503L else 200L, content = calls)
  })
  first <- cache_static_app_ui(original)
  second <- cache_static_app_ui(original)
  req <- list(REQUEST_METHOD = "GET", PATH_INFO = "/", QUERY_STRING = "")
  expect_null(first$httpHandler(req))
  expect_equal(first$httpHandler(req)$status, 503L)
  expect_equal(first$httpHandler(req)$content, 3L)
  expect_equal(first$httpHandler(req)$content, 3L)
  expect_equal(second$httpHandler(req)$content, 4L)
})
