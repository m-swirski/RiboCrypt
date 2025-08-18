sample_info_ui <- function(id, label = "sample_info") {
  ns <- NS(id)
  tabPanel(title = "Samples", icon = icon("rectangle-list"),
           h2("All samples (processed subset)"),
           mainPanel(
             DT::DTOutput(ns("sample_info")) #%>% shinycssloaders::withSpinner(color="#0dc5c1")
           )
  )
}

sample_info_server <- function(id, metadata) {
  moduleServer(
    id,
    function(input, output, session) {
      tableData <- reactiveVal(metadata)
      output$sample_info <- DT::renderDT(tableData(),
                                         extensions = 'Buttons',
                                         filter = "top",
                                         options = list(dom = 'Bfrtip',
                                                        buttons = NULL),
                                         server = TRUE)
      
      proxyTable <- DT::dataTableProxy("sample_info")
      
    }
  )
}
