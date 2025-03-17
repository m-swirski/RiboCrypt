run_info_ui <- function(id, label = "run_info") {
  ns <- NS(id)
  tabPanel(title = "Runs", icon = icon("rectangle-list"),
           h2("All runs (processed & processing failed)"),
           mainPanel(
             DT::DTOutput(ns("run_info")) %>% shinycssloaders::withSpinner(color="#0dc5c1")
           )
  )
}

run_info_server <- function(id, metadata) {
  moduleServer(
    id,
    function(input, output, session) {
      columns_to_show <- c("study_accession", "Run", "ScientificName", "sample_title", "LIBRARYTYPE", "REPLICATE", "CONDITION", "INHIBITOR",
                           "BATCH", "TIMEPOINT", "TISSUE", "CELL_LINE", "GENE", "FRACTION")
      metadata_subset <- metadata[, colnames(metadata) %in% columns_to_show, with = FALSE]
      output$run_info <- DT::renderDT(metadata_subset,
                                      extensions = 'Buttons',
                                      filter = "top",
                                      options = list(dom = 'Bfrtip',
                                                     buttons = NULL))
    }
  )
}
