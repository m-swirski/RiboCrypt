sra_search_ui <- function(id, label = "sra_search") {
  ns <- NS(id)
  tabPanel(title = "SRA search", icon = icon("rectangle-list"),
           h2("Metadata search page"),
           sidebarLayout(
             sidebarPanel(
               textInput(ns("accession"), "Study accession number (SRP/GEO/PRJNA)",
                         "GSE13750"),
               actionButton(ns("go"), "Search", icon = icon("magnifying-glass")),
             ),
             mainPanel(
               textOutput(ns("abstract")),
               dataTableOutput(ns("metadata")) %>% shinycssloaders::withSpinner(color="#0dc5c1")
             )
           )
  )
}

sra_search_server <- function(id) {
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
