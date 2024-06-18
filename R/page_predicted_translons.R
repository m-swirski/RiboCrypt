predicted_translons_ui <- function(id, all_exp_translons, label = "predicted_translons") {
  ns <- NS(id)
  tabPanel(title = "Predicted Translons", icon = icon("rectangle-list"),
           h2("Predicted Translons overview"),
           sidebarLayout(
             sidebarPanel(
               experiment_input_select(all_exp_translons$name, ns),
               actionButton(ns("go"), "Search", icon = icon("magnifying-glass"))
             ),
             mainPanel(
               DT::DTOutput(ns("translon_table")) %>% shinycssloaders::withSpinner(color="#0dc5c1")
             )
           )
  )
}

predicted_translons_server <- function(id) {
  moduleServer(
    id,
    function(input, output, session) {
      md <- eventReactive(input$go,{
        species <- isolate(input$dff)
        df <- read.experiment(species, validate = FALSE)
        table_path <- file.path(dirname(df@fafile),
                                "predicted_translons",
                                "predicted_translons_with_sequence.fst")
        if (file.exists(table_path)) {
          translon_table <- fst::read_fst(table_path)
        }

        reactiveValues(translon_table = translon_table)
      })
      output$translon_table <- DT::renderDT(md()$translon_table,
                                            extensions = 'Buttons',
                                            filter = "top",
                                            options = list(dom = 'Bfrtip',
                                                           buttons = c('csv', 'excel'))
      )
    }
  )
}
