test_that("Predicted Translons keeps its experiment filter and module IDs", {
  experiments <- data.table::data.table(
    name = c("human_all_merged_l50", "mouse_all_merged", "Escherichia_coli_all_merged",
             "human_study", "human_all_merged_rna"),
    libtypes = c("RFP", "RFP", "RFP", "RFP", "RNA"), organism = "test")
  selected <- predicted_translons_experiments(experiments)
  expect_equal(selected$name, c("human_all_merged_l50", "mouse_all_merged"))
  tab <- predicted_translons_ui("predicted_translons", selected)
  expect_identical(tab$attribs$title, "Translons")
  html <- htmltools::renderTags(tab)$html
  expect_match(html, 'data-value="Predicted Translons"', fixed = TRUE)
  expect_match(html, 'id="predicted_translons-go"', fixed = TRUE)
  expect_match(html, 'id="predicted_translons-translon_table"', fixed = TRUE)
  metadata <- htmltools::renderTags(metadata_ui("metadata", experiments, experiments[0]))$html
  expect_false(grepl("predicted_translons", metadata, fixed = TRUE))
})

test_that("Metadata server no longer registers the moved translon module", {
  calls <- character()
  record <- function(id, ...) calls <<- c(calls, id)
  testthat::local_mocked_bindings(sample_info_server = record, study_info_server = record,
    sra_search_server = record, umap_server = record,
    predicted_translons_server = function(...) stop("Must initialize independently"),
    .package = "RiboCrypt")
  metadata_server("metadata", data.table::data.table(), data.frame(), data.table::data.table(),
                  c(search_on_init = "FALSE"))
  expect_equal(calls, c("sample_info", "study_info", "sra_search", "umap"))
})
