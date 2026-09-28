test_that("simplified columns preserve display order and optional localization columns", {
  columns <- c("link", "ID", "ensembl_gene_id", "alignment", "ensembl_tx_name",
               "external_gene_name", "type", "length", "sequence_AA", "mapability_score",
               "biotype", "LOC_nucleus", "LOC_mitochondria", "other_LOC_value")
  expect_equal(columns[translon_simplified_columns(columns)],
    c("link", "ID", "ensembl_gene_id", "external_gene_name", "type", "length",
      "mapability_score", "biotype", "LOC_nucleus", "LOC_mitochondria"))
  expect_equal(translon_simplified_columns(c("a", "b")), 1:2)
  expect_length(translon_simplified_columns(character()), 0)
  expect_equal(translon_simplified_columns(c("a", "b", "c", "unwanted")), 1:3)
})

test_that("table controls use post-link column indices and preserve ID clicks and exports", {
  data <- data.table::data.table(ensembl_gene_id = "GENE", alignment = "1:1-9:+",
    ensembl_tx_name = "TX", external_gene_name = "G", type = "uORF", length = 9,
    ID = "pep1", coordinates = "1:9", mapability_score = 1, biotype = "protein_coding",
    LOC_nucleus = 0.8)
  data.table::setattr(data, "exp", "experiment")
  before <- data.table::copy(data)
  session <- list(ns = shiny::NS("predicted_translons"), clientData = list(
    url_hostname = "localhost", url_port = "7821", url_pathname = "/", url_protocol = "http:"))
  widget <- render_translon_datatable(data, session)
  expect_identical(names(widget$x$data)[1:3], c("link", "ID", "ensembl_gene_id"))
  expect_match(widget$x$data$link, "GENE", fixed = TRUE)
  expect_true(widget$x$options$scrollX)
  expect_match(as.character(widget$x$callback), "predicted_translons-simplified", fixed = TRUE)
  expect_match(as.character(widget$x$callback), "predicted_translons-translon_id_click", fixed = TRUE)
  for (button in widget$x$options$buttons) expect_equal(button$exportOptions$columns, ":visible")
  expect_identical(data, before)
})

test_that("the toolbar includes all controls above the full-width table", {
  html <- htmltools::renderTags(predicted_translons_ui("p", data.table::data.table(name = "exp")))$html
  expect_match(html, 'class="rc-translon-toolbar"', fixed = TRUE)
  for (group in c("study-group", "simplified", "downloads")) {
    expect_match(html, paste0('class="rc-translon-', group, '"'), fixed = TRUE)
  }
  expect_match(html, "rc-translon-excel", fixed = TRUE)
  for (id in c("dff", "go", "simplified", "trigger_download_csv", "trigger_download_excel")) {
    expect_lt(regexpr(paste0('id="p-', id, '"'), html, fixed = TRUE)[1],
              regexpr('id="p-translon_table"', html, fixed = TRUE)[1])
  }
  expect_false(grepl('class="well"', html, fixed = TRUE))
  expect_match(html, 'id="p-simplified" type="checkbox"', fixed = TRUE)
})

test_that("the shipped column toggle preserves table state and rebinds safely", {
  node <- Sys.which("node")
  skip_if(!nzchar(node), "Node.js is required for JavaScript callback tests")
  withr::local_envvar(RIBOCRYPT_TRANSLON_JS = system.file("js", "translon_table_controls.js", package = "RiboCrypt"))
  result <- suppressWarnings(system2(node, c("--test", shQuote(test_path("js", "translon-table-controls.cjs"))),
                                     stdout = TRUE, stderr = TRUE))
  expect_null(attr(result, "status"), info = paste(result, collapse = "\n"))
})
