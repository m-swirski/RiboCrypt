sample_info_ui <- function(id, label = "sample_info") {
  ns <- NS(id)
  tabPanel(title = "Samples", icon = icon("rectangle-list"),
           h2("All samples (processed subset)"),
           mainPanel(
             DT::DTOutput(ns("sample_info")) %>% shinycssloaders::withSpinner(color="#0dc5c1")
           )
  )
}

sample_info_server <- function(id, metadata, search_on_init = "") {
  moduleServer(
    id,
    function(input, output, session) {
      colnames(metadata)[colnames(metadata) == "Run"] <- "Sample"
      colnames(metadata)[colnames(metadata) == "ScientificName"] <- "Organism"
      output$sample_info <- DT::renderDT(metadata,
                                         extensions = 'Buttons',
                                         filter = "top",
                                         options = list(dom = 'Bfrtip',
                                                        buttons = NULL,
                                                        search = list(search = as.character(search_on_init))))
    }
  )
}
