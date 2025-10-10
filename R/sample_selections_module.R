sampleSelectionsUi <- function(id) {
  ns <- NS(id)
  fluidRow(
    fluidRow(
      column(
        2,
        selectizeInput(ns("activeSelectionSelect"), "Selection", choices = list("New"))
      )
    ),
    fluidRow(
      column(
        12,
        uiOutput(ns("sampleSelectionsUi"))
      )
    )
  )
}

newSelectionUi <- function(ns) {
  fluidRow(
    column(2, actionButton(ns("saveAsSelection"), "Save selection"))
  )
}

selectionUi <- function(ns) {
  fluidRow(
    sampleTableUi(ns("activeSelectionTable"))
  )
}

sampleSelectionsServer <- function(id, metadata, rInitialSelection) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # reactive values
    counter <- reactiveVal(1)
    rSelections <- reactiveVal(list())
    rActiveSelectionId <- reactiveVal()
    rActiveSelection <- reactive({
      x <- base::Find(function(elem) {
        elem$id == rActiveSelectionId()
      }, rSelections())
      x$rSelection
    }) %>% bindEvent(rActiveSelectionId())
    selectionServer <- sampleTableServer("activeSelectionTable", metadata, rActiveSelection)

    # Observer for handling adding a new selection
    observe({
      count <- counter()
      counter(count + 1)

      selections <- rSelections()
      rSelection <- reactiveVal(rInitialSelection())
      newSelection <- list(
        id = as.character(count),
        label = as.character(count),
        rSelection = rSelection
      )

      rSelections(
        c(selections, list(newSelection))
      )
      rActiveSelectionId(newSelection$id)
    }) %>% bindEvent(input$saveAsSelection)

    # Observers for handling interaction with the select input
    observe({
      updateSelectizeInput(
        session,
        "activeSelectionSelect",
        choices = c(lapply(rSelections(), function(elem) {
          elem$id
        }), list("New")),
        selected = rActiveSelectionId()
      )
    }) %>% bindEvent(rSelections())

    observe({
      rActiveSelectionId(input$activeSelectionSelect)
    }) %>% bindEvent(input$activeSelectionSelect)

    # Observers for handling interactions with the outside world
    observe({
      req(!is.null(rActiveSelection()))
      rActiveSelection()(rInitialSelection())
    }) %>% bindEvent(rInitialSelection())

    observe({
      req(!is.null(rActiveSelection()))
      message <- {
        curveIndex <- rActiveSelection()()$curveIndex
        pointIndex <- rActiveSelection()()$pointIndex

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
    }) %>% bindEvent(rActiveSelection())

    # Dynamic Ui
    output$sampleSelectionsUi <- renderUI({
      if (is.null(rActiveSelection())) {
        newSelectionUi(ns)
      } else {
        selectionUi(ns)
      }
    })

    list(
      selections = rSelections,
      activeSelectionId = rActiveSelectionId,
      activeSelection = rActiveSelection
    )
  })
}
