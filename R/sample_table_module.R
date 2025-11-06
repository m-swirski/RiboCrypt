sampleTableUi <- function(id) {
  ns <- shiny::NS(id)
  DT::DTOutput(ns("sampleTable"))
}

sampleTableServer <- function(id, metadata, rSelection, rFilteredSelection) {
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

    shiny::observe({
      selectedInTable <- tableData()[input$sampleTable_rows_selected]$Sample
      if (is.null(selectedInTable)) {
        rFilteredSelection(NULL)
      } else {
        filteredSelection <- list(
          sample = rSelection()$sample[rSelection()$sample %in% selectedInTable],
          curveIndex = rSelection()$curveIndex[rSelection()$sample %in% selectedInTable],
          pointIndex = rSelection()$pointIndex[rSelection()$sample %in% selectedInTable]
        )
        rFilteredSelection(filteredSelection)
      }
    }) %>% shiny::bindEvent(rSelection(), input$sampleTable_rows_selected)
  })
}
