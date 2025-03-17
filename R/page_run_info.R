sample_info_ui <- function(id, label = "sample_info") {
  ns <- NS(id)
  tabPanel(title = "Samples", icon = icon("rectangle-list"),
           h2("All samples (processed subset)"),
           mainPanel(
             DT::DTOutput(ns("sample_info")) %>% shinycssloaders::withSpinner(color="#0dc5c1")
           )
  )
}

sample_info_server <- function(id, metadata) {
  moduleServer(
    id,
    function(input, output, session) {
      columns_to_show <- c("study_accession", "Run", "ScientificName", "sample_title", "BioProject",
                           "LIBRARYTYPE", "REPLICATE", "CONDITION", "INHIBITOR",
                           "BATCH", "TIMEPOINT", "TISSUE", "CELL_LINE", "GENE", "FRACTION")
      metadata_subset <- metadata[, colnames(metadata) %in% columns_to_show, with = FALSE]
      colnames(metadata_subset)[colnames(metadata_subset) == "Run"] <- "Sample"

      output$sample_info <- DT::renderDT(metadata_subset,
                                         extensions = 'Buttons',
                                         filter = "top",
                                         options = list(dom = 'Bfrtip',
                                                        buttons = NULL))
    }
  )
}
