test_that("input_to_list drops ignored inputs and adds user info", {
  input <- shiny::reactiveValues(
    gene = "GENE1",
    tx = "TX1",
    go = 1,
    toggle_settings = 0,
    select_all_btn = 0,
    c__shinyjquiBookmarkState__resizable = NULL,
    c_is_resizing = FALSE,
    c_size = NULL
  )
  res <- shiny::isolate(RiboCrypt:::input_to_list(
    input,
    list(id = "user-1", browser = "firefox")
  ))
  expect_true("gene" %in% names(res))
  expect_true("tx" %in% names(res))
  expect_false("go" %in% names(res))
  expect_false("toggle_settings" %in% names(res))
  expect_false("select_all_btn" %in% names(res))
  expect_true("user.info.browser" %in% names(res))
  expect_equal(res[["user.info.browser"]], "firefox")
})

test_that("go_when_input_is_ready triggers kickoff when inputs match", {
  browser_options <- c(
    plot_on_start = "TRUE",
    default_gene = "GENE1",
    default_isoform = "TX1",
    default_libs = "A,B"
  )
  input <- list(
    gene = "GENE1",
    tx = "TX1",
    library = c("A", "B")
  )
  libs <- shiny::reactiveVal(c("A", "B"))
  fired <- shiny::reactiveVal(FALSE)
  kickoff <- shiny::reactiveVal(FALSE)

  shiny::isolate(RiboCrypt:::go_when_input_is_ready(input, browser_options, fired, kickoff, libs))
  expect_true(shiny::isolate(isTRUE(fired())))
  expect_true(shiny::isolate(isTRUE(kickoff())))
})

test_that("go_when_input_is_ready does not trigger when plot_on_start is FALSE", {
  browser_options <- c(
    plot_on_start = "FALSE",
    default_gene = "GENE1",
    default_isoform = "TX1",
    default_libs = "A,B"
  )
  input <- list(
    gene = "GENE1",
    tx = "TX1",
    library = c("A", "B")
  )
  libs <- shiny::reactiveVal(c("A", "B"))
  fired <- shiny::reactiveVal(FALSE)
  kickoff <- shiny::reactiveVal(FALSE)

  shiny::isolate(RiboCrypt:::go_when_input_is_ready(input, browser_options, fired, kickoff, libs))
  expect_true(shiny::isolate(isTRUE(fired())))
  expect_false(shiny::isolate(isTRUE(kickoff())))
})

test_that("browser_ui returns a shiny tag object", {
  all_exp <- data.frame(
    organism = c("org1", "org2"),
    name = c("exp1", "exp2"),
    stringsAsFactors = FALSE
  )
  browser_options <- c(
    default_gene = "GENE1",
    default_isoform = "TX1",
    default_libs = "A|B",
    default_view_mode = "tx",
    collapsed_introns_width = "100",
    full_annotation = "FALSE",
    translons = "FALSE",
    translons_transcode = "FALSE",
    hide_settings = "TRUE",
    default_frame_type = "all",
    default_kmer = "1"
  )
  gene_names_init <- data.frame(label = "GENE1", stringsAsFactors = FALSE)
  libs <- c("A", "B")

  ui <- shiny::isolate(RiboCrypt:::browser_ui("test", all_exp, browser_options, gene_names_init, libs))
  expect_true(inherits(ui, "shiny.tag") || inherits(ui, "shiny.tag.list"))
})
