createNewSelectionChoice <- "New selection..."

sampleSelectionsUi <- function(id) {
  ns <- shiny::NS(id)
  shiny::fluidRow(
    shiny::column(
      2,
      shiny::selectizeInput(ns("activeSelectionSelect"), "Selection", choices = list())
    )
  )
}

sampleSelectionsServer <- function(id, metadata, rSelection, rFilteredSelection) {
  shiny::moduleServer(id, function(input, output, session) {
    # reactive values
    counter <- shiny::reactiveVal(2)
    rSelections <- {
      selections <- list()
      selections[[as.character(1)]] <- NULL
      filteredSelections <- list()
      filteredSelections[[as.character(1)]] <- NULL
      shiny::reactiveVal(list(
        index = c(as.character(1)),
        selections = selections,
        filteredSelections = filteredSelections
      ))
    }

    rActiveSelectionId <- shiny::reactiveVal(as.character(1))

    rActiveSelection <- shiny::reactive({
      shiny::req(!is.null(rActiveSelectionId()) && rActiveSelectionId() != "")
      rSelections()$selections[[rActiveSelectionId()]]
    }) %>% shiny::bindEvent(rSelections(), rActiveSelectionId())

    rActiveFilteredSelection <- shiny::reactive({
      shiny::req(!is.null(rActiveSelectionId()) && rActiveSelectionId() != "")
      rSelections()$filteredSelections[[rActiveSelectionId()]]
    }) %>% shiny::bindEvent(rSelections(), rActiveSelectionId())

    # Observer for creating a new selection
    shiny::observe({
      shiny::req(input$activeSelectionSelect == createNewSelectionChoice)
      newSelectionId <- counter()
      counter(newSelectionId + 1)

      selections <- rSelections()

      selections$index <- c(selections$index, as.character(newSelectionId))
      selections$selections[[as.character(newSelectionId)]] <- NULL
      selections$filteredSelections[[as.character(newSelectionId)]] <- NULL

      rSelections(selections)
      rActiveSelectionId(as.character(newSelectionId))
    }) %>% shiny::bindEvent(input$activeSelectionSelect)

    # Observers for selecting the active selection
    shiny::observe({
      shiny::updateSelectizeInput(
        session,
        "activeSelectionSelect",
        choices = c(rSelections()$index, createNewSelectionChoice),
        selected = rActiveSelectionId()
      )
    }) %>% shiny::bindEvent(rSelections())

    shiny::observe({
      shiny::req(
        !is.null(input$activeSelectionSelect) &&
          input$activeSelectionSelect != "" &&
          input$activeSelectionSelect != createNewSelectionChoice
      )
      rActiveSelectionId(input$activeSelectionSelect)
    }) %>% shiny::bindEvent(input$activeSelectionSelect)

    # Observers for handling interactions with selection
    shiny::observe({
      selections <- rSelections()
      selections$selections[[rActiveSelectionId()]] <- rSelection()

      rSelections(selections)
    }) %>% shiny::bindEvent(rSelection())

    shiny::observe({
      rSelection(rActiveSelection())
    }) %>% shiny::bindEvent(rActiveSelection(), ignoreNULL = FALSE)

    shiny::observe({
      message <- if (is.null(rActiveSelection())) {
        list(
          curveIndex = list(),
          pointIndex = list()
        )
      } else {
        curveIndex <- rActiveSelection()$curveIndex
        pointIndex <- rActiveSelection()$pointIndex

        list(
          curveIndex = unique(curveIndex),
          pointIndex = lapply(unique(curveIndex), function(idx) {
            pointIndex[curveIndex == idx]
          })
        )
      }

      session$sendCustomMessage(
        "samplesActiveSelectionChanged",
        message
      )
    }) %>% shiny::bindEvent(rActiveSelection(), ignoreNULL = FALSE)

    # Observers for handling interactions with filtered selection

    shiny::observe({
      selections <- rSelections()

      if (is.null(rFilteredSelection())) {
        selections$filteredSelections[[rActiveSelectionId()]] <- rSelection()
      } else {
        selections$filteredSelections[[rActiveSelectionId()]] <- rFilteredSelection()
      }

      rSelections(selections)
    }) %>% shiny::bindEvent(rSelection(), rFilteredSelection(), ignoreNULL = FALSE)

    shiny::observe({
      shiny::req(!is.null(rActiveFilteredSelection()))
      message <- {
        curveIndex <- rActiveFilteredSelection()$curveIndex
        pointIndex <- rActiveFilteredSelection()$pointIndex

        list(
          curveIndex = unique(curveIndex),
          pointIndex = lapply(unique(curveIndex), function(idx) {
            pointIndex[curveIndex == idx]
          })
        )
      }

      session$sendCustomMessage(
        "samplesActiveFilteredSelectionChanged",
        message
      )
    }) %>% shiny::bindEvent(rActiveFilteredSelection())

    rSelections
  })
}
