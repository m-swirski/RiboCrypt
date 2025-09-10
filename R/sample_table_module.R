sampleTableUi <- function(id) {
  ns <- NS(id)
  fluidRow(uiOutput(ns("sampleTableUi")))
}


sampleTableServer <- function(id, metadata, rInitialSelectedSamples) {
  moduleServer(id, function(input, output, session){
    ns <- session$ns
    selectedSamples <- reactiveVal()
    tableData <- reactive({
      req(!is.null(selectedSamples()))
      metadata[Sample %in% selectedSamples()]
    })
    
    observe({
      req(rInitialSelectedSamples)
      req(is.null(selectedSamples()))
      selectedSamples(rInitialSelectedSamples)
    }) %>% bindEvent(input$saveAsGroup)
    
    observe({
      req(!is.null(selectedSamples()))
      selectedSamples(NULL)
    }) %>% bindEvent(input$clearSelection)
    
    observe({
      req(!is.null(selectedSamples()))
      indexesToRemove <- input$sampleTable_rows_selected
      samplesToRemove <- tableData()[indexesToRemove][["Sample"]]
      currentSelection <- selectedSamples()
      selectedSamples(currentSelection[!currentSelection %in% samplesToRemove])
    }) %>% bindEvent(input$removeSelected)
    
    output$sampleTable <-
      DT::renderDT(tableData(), filter = "top", options = list(dom = 'Bfrtip'))
    
    output$sampleTableUi <- renderUI({
      if (is.null(selectedSamples())) {
        actionButton(ns("saveAsGroup"), "Save selection")
      } else {
        fluidRow(
          column(10, DT::DTOutput(ns("sampleTable"))),
          column(2, 
                 actionButton(ns("clearSelection"), "Clear"),
                 actionButton(ns("removeSelected"), "Remove"))
          )
      }
    })
    
    return(selectedSamples)
  })
}
