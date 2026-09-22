test_that("hidden servers initialize once and retain their reactive state", {
  starts <- 0L
  shiny::testServer(function(input, output, session) {
    on_first_tab(input, c("MegaBrowser", "Observatory"), function() {
      starts <<- starts + 1L
      shiny::moduleServer("deferred", function(input, output, session) {
        value <- shiny::reactiveVal(0L)
        shiny::observeEvent(input$increment, value(value() + 1L))
        output$value <- shiny::renderText(value())
      })
    })
  }, {
    session$setInputs(navbarID = "browser")
    expect_equal(starts, 0L)
    session$setInputs(navbarID = "Observatory")
    expect_equal(starts, 1L)
    session$setInputs(`deferred-increment` = 1)
    expect_equal(output[["deferred-value"]], "1")
    session$setInputs(navbarID = "browser")
    session$setInputs(navbarID = "MegaBrowser", `deferred-increment` = 2)
    expect_equal(output[["deferred-value"]], "2")
    expect_equal(starts, 1L)
  })
})

test_that("initial tab selection initializes the server in each new session", {
  starts <- 0L
  server <- function(input, output, session) {
    on_first_tab(input, "Studies", function() starts <<- starts + 1L)
  }
  for (i in 1:2) {
    shiny::testServer(server, {
      session$setInputs(navbarID = "Studies")
      expect_equal(starts, i)
      session$setInputs(navbarID = "browser")
      session$setInputs(navbarID = "Studies")
      expect_equal(starts, i)
    })
  }
})

test_that("gene validation preserves selection with repeated transcript labels", {
  genes <- data.table::data.table(
    value = paste0("TX", 1:5), label = c("B", "B", NA, "A", "A")
  )
  expect_true(gene_exists_in_gene_list(genes, "A"))
  expect_false(gene_exists_in_gene_list(genes, NA_character_))
  expect_false(gene_exists_in_gene_list(genes, "absent"))
  expect_identical(resolve_gene_selection(genes, "A", "B"), "A")
  expect_identical(resolve_gene_selection(genes, "absent", "A"), "A")
  expect_identical(resolve_gene_selection(genes, "absent"), "B")
  expect_identical(resolve_gene_selection(genes[0]), character())
})
