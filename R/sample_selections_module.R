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

sampleSelectionsServer <- function(id, metadata, rPrimarySelection, rSecondarySelection) {
  shiny::moduleServer(id, function(input, output, session) {
    # reactive values
    rSelections <- {
      selections <- list()
      selections[[as.character(1)]] <- NULL
      shiny::reactiveVal(list(
        index = c(as.character(1)),
        selections = selections
      ))
    }

    rActiveSelectionId <- shiny::reactiveVal(as.character(1))
    rActiveSelection <- shiny::reactive({
      shiny::req(!is.null(rActiveSelectionId()) && rActiveSelectionId() != "")
      rSelections()$selections[[rActiveSelectionId()]]
    }) %>% shiny::bindEvent(rSelections(), rActiveSelectionId())
    counter <- shiny::reactiveVal(2)

    # Observers for handling interaction with the select input
    shiny::observe({
      shiny::updateSelectizeInput(
        session,
        "activeSelectionSelect",
        choices = c(rSelections()$index, createNewSelectionChoice),
        selected = rActiveSelectionId()
      )
    }) %>% shiny::bindEvent(rSelections())

    shiny::observe({
      shiny::req(!is.null(input$activeSelectionSelect) && input$activeSelectionSelect != "")
      rActiveSelectionId(input$activeSelectionSelect)
    }) %>% shiny::bindEvent(input$activeSelectionSelect)

    # Observer for handling creating a new selection
    shiny::observe({
      shiny::req(input$activeSelectionSelect == createNewSelectionChoice)
      newSelectionId <- counter()
      counter(newSelectionId + 1)

      selections <- rSelections()

      selections$index <- c(selections$index, as.character(newSelectionId))
      selections$selections[[as.character(newSelectionId)]] <- NULL

      rSelections(selections)
      rActiveSelectionId(as.character(newSelectionId))
    }) %>% shiny::bindEvent(input$activeSelectionSelect)

    # Observers for handling interactions with the outside world
    shiny::observe({
      selections <- rSelections()
      selections$selections[[rActiveSelectionId()]] <- rPrimarySelection()

      rSelections(selections)
    }) %>% shiny::bindEvent(rPrimarySelection())

    shiny::observe({
      shiny::req(!is.null(rActiveSelectionId()) && rActiveSelectionId() != "")
      rSecondarySelection(rActiveSelection())
    }) %>% shiny::bindEvent(rActiveSelection(), ignoreNULL = FALSE)

    # observer({
    #   req(!is.null(rActiveSelection()))
    #   rActiveSelection()(rSecondarySelection())
    # }) %>% bindEvent(rSecondarySelection())

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
        "samplesSelectionChanged",
        message
      )
    }) %>% shiny::bindEvent(rActiveSelection(), ignoreNULL = FALSE)

    rSelections
  })
}
