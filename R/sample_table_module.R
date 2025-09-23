sampleTableUi <- function(id) {
  ns <- NS(id)
  DT::DTOutput(ns("sampleTable"))
}


sampleTableServer <- function(id, metadata, initialSelectedSamples) {
  moduleServer(id, function(input, output, session){
    ns <- session$ns
    selectedSamples <- reactiveVal(initialSelectedSamples)
    
    tableData <- reactive({
      req(!is.null(selectedSamples()))
      metadata[Sample %in% selectedSamples()]
    })
    
    
    removeObserver <- observe({
      req(!is.null(selectedSamples()))
      indexesToRemove <- input$sampleTable_rows_selected
      samplesToRemove <- tableData()[indexesToRemove][["Sample"]]
      currentSelection <- selectedSamples()
      selectedSamples(currentSelection[!currentSelection %in% samplesToRemove])
    }) %>% bindEvent(input$removeSelected)
    
    output$sampleTable <-
      DT::renderDT(tableData(), filter = "top", options = list(dom = 'Bfrtip'))
    
    return(list(id = id, observers = c(removeObserver), selectedSamples = selectedSamples))
  })
}
