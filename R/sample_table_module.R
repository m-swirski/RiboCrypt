sampleTableUi <- function(id) {
  ns <- NS(id)
  DT::DTOutput(ns("sampleTable"))
}


sampleTableServer <- function(id, metadata, rSelection) {
  moduleServer(id, function(input, output, session){
    ns <- session$ns
    selectedSamples <- reactive({
      req(!is.null(rSelection()))
      rSelection()()$sample
    })
    
    tableData <- reactive({
      req(!is.null(selectedSamples()))
      metadata[Sample %in% selectedSamples()]
    })
    
    observe({
      req(!is.null(selectedSamples()))
      indexesToRemove <- input$sampleTable_rows_selected
      samplesToRemove <- tableData()[indexesToRemove][["Sample"]]
      currentSelection <- rSelection()()
      rSelection()(currentSelection[!currentSelection$sample %in% samplesToRemove])
    }) %>% bindEvent(input$removeSelected)
    
    output$sampleTable <-
      DT::renderDT(tableData(), filter = "top", options = list(dom = 'Bfrtip'))
  })
}
