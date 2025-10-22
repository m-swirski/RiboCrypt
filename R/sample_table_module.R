sampleTableUi <- function(id) {
  ns <- NS(id)
  DT::DTOutput(ns("sampleTable"))
}

sampleTableServer <- function(id, metadata, rSelection) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    selectedSamples <- reactive({
      req(!is.null(rSelection()))
      rSelection()$sample
    })

    tableData <- reactive({
      req(!is.null(selectedSamples()))
      metadata[Sample %in% selectedSamples()]
    })

    output$sampleTable <-
      DT::renderDT(tableData(), filter = "top", options = list(dom = "Bfrtip"))
  })
}
