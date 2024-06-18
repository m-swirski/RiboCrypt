study_info_ui <- function(id, label = "study_info") {
  ns <- NS(id)
  tabPanel(title = "Studies", icon = icon("rectangle-list"),
           h2("Supported studies overview"),
           mainPanel(
             DT::DTOutput(ns("metadata")) %>% shinycssloaders::withSpinner(color="#0dc5c1")
           )
  )
}

study_info_server <- function(id, all_exp) {
  moduleServer(
    id,
    function(input, output, session) {
      output$metadata <- DT::renderDT(all_exp,
                                      extensions = 'Buttons',
                                      filter = "top",
                                      options = list(dom = 'Bfrtip',
                                                     buttons = c('csv', 'excel')))
    }
  )
}
