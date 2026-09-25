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

test_that("collection tabs initialize independently and retain their shared state", {
  for (first in c("Observatory", "MegaBrowser")) {
    starts <- character()
    shiny::testServer(function(input, output, session) {
      on_first_tab(input, c("MegaBrowser", "Observatory"), function() {
        starts <<- c(starts, "shared")
        selection <- shiny::reactiveVal("initial")
        on_first_tab(input, "MegaBrowser", function() {
          starts <<- c(starts, "MegaBrowser")
          shiny::observeEvent(input$mega_edit, selection(input$mega_edit))
          output$mega <- shiny::renderText(selection())
        })
        on_first_tab(input, "Observatory", function() {
          starts <<- c(starts, "Observatory")
          shiny::observeEvent(input$obs_edit, selection(input$obs_edit))
          output$obs <- shiny::renderText(selection())
        })
      })
    }, {
      session$setInputs(navbarID = "browser")
      expect_length(starts, 0L)
      session$setInputs(navbarID = first)
      expect_identical(starts, c("shared", first))
      other <- setdiff(c("Observatory", "MegaBrowser"), first)
      session$setInputs(navbarID = other)
      expect_identical(starts, c("shared", first, other))
      session$setInputs(obs_edit = "selected")
      expect_identical(output$mega, "selected")
      session$setInputs(navbarID = "browser")
      session$setInputs(navbarID = first, mega_edit = "changed")
      expect_identical(output$obs, "changed")
      expect_identical(starts, c("shared", first, other))
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
