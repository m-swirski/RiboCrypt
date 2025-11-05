sampleTableUi <- function(id) {
  ns <- shiny::NS(id)
  DT::DTOutput(ns("sampleTable"))
}

sampleTableServer <- function(id, metadata, rSelection) {
  shiny::moduleServer(id, function(input, output, session) {
    selectedSamples <- shiny::reactive({
      if (is.null(rSelection())) {
        c("")
      } else {
        rSelection()$sample
      }
    })

    tableData <- shiny::reactive({
      shiny::req(!is.null(selectedSamples()))
      metadata[Sample %in% selectedSamples()]
    })

    output$sampleTable <-
      DT::renderDT(tableData(), filter = "top", options = list(dom = "Bfrtip"))
  })
}
