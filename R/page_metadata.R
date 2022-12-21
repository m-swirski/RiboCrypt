metadata_ui <- function(id, label = "metadata") {
  ns <- NS(id)
  tabPanel(title = "metadata", icon = icon("rectangle-list"),
           h2("Metadata search page"),
           sidebarLayout(
             sidebarPanel(
               textInput(ns("accession"), "Study accession number (SRP/GEO/PRJNA)"),
               actionButton(ns("go"), "Plot", icon = icon("rocket")),
             ),
             mainPanel(
               textOutput(ns("abstract")),
               dataTableOutput(ns("metadata"))
             )
           )
  )
}

metadata_server <- function(id) {
  moduleServer(
    id,
    function(input, output, session) {
      md <- eventReactive(input$go,{
        abstract <- capture.output({
          meta_dt <- download.SRA.metadata(input$accession)
          # TODO: Make a more clever shower
          meta_dt <- setcolorder(meta_dt, rev(seq(ncol(dt))))
        })
        reactiveValues(abstract = abstract,
                       meta_dt = meta_dt)
      })
      output$abstract <- renderText(md()$abstract)
      output$metadata <- renderDataTable(md()$meta_dt)
    }
  )
}
