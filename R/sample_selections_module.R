createNewSelectionChoice <- "New selection..."

sampleSelectionsUi <- function(id) {
  ns <- NS(id)
  fluidRow(
    column(
      2,
      selectizeInput(ns("activeSelectionSelect"), "Selection", choices = list(createNewSelectionChoice))
    )
  )
}

sampleSelectionsServer <- function(id, metadata, rPrimarySelection, rSecondarySelection) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # reactive values
    counter <- reactiveVal(2)
    rSelections <- reactiveVal(
      list(
        list(
          id = as.character(1),
          label = as.character(1),
          rSelection = reactiveVal(NULL)
        )
      )
    )
    rActiveSelectionId <- reactiveVal(as.character(1))
    rActiveSelection <- reactive({
      x <- base::Find(function(elem) {
        elem$id == rActiveSelectionId()
      }, rSelections())
      x$rSelection()
    }) %>% bindEvent(rActiveSelectionId())

    # Observers for handling interaction with the select input
    observe({
      updateSelectizeInput(
        session,
        "activeSelectionSelect",
        choices = c(lapply(rSelections(), function(elem) {
          elem$id
        }), list(createNewSelectionChoice)),
        selected = rActiveSelectionId()
      )
    }) %>% bindEvent(rSelections())

    observe({
      rActiveSelectionId(input$activeSelectionSelect)
    }) %>% bindEvent(input$activeSelectionSelect)

    # Observer for handling creating a new selection
    observe({
      req(input$activeSelectionSelect == createNewSelectionChoice)
      count <- counter()
      counter(count + 1)

      selections <- rSelections()
      rSelection <- reactiveVal(rPrimarySelection())
      newSelection <- list(
        id = as.character(count),
        label = as.character(count),
        rSelection = rSelection
      )

      rSelections(
        c(selections, list(newSelection))
      )
      rActiveSelectionId(newSelection$id)
    }) %>% bindEvent(input$activeSelectionSelect)

    # Observers for handling interactions with the outside world
    observe({
      req(!is.null(rActiveSelection()))
      rSecondarySelection(rActiveSelection())
    }) %>% bindEvent(rActiveSelection())

    observe({
      req(!is.null(rActiveSelection()))
      rActiveSelection(rPrimarySelection())
    }) %>% bindEvent(rPrimarySelection())

    # observer({
    #   req(!is.null(rActiveSelection()))
    #   rActiveSelection()(rSecondarySelection())
    # }) %>% bindEvent(rSecondarySelection())

    observe({
      req(!is.null(rActiveSelection()))
      message <- {
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
    }) %>% bindEvent(rActiveSelection())

    list(
      selections = rSelections,
      activeSelectionId = rActiveSelectionId,
      activeSelection = rActiveSelection
    )
  })
}
