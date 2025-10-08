sampleSelectionsUi <- function(id) {
  ns <- NS(id)
  fluidRow(
    column(
      10,
      uiOutput(ns("sampleSelectionsUi"))
    ),
    column(
      2,
      selectizeInput(ns("activeSelectionSelect"), "Selection", choices = list("New"))
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
  moduleServer(id, function(input, output, session){
    ns <- session$ns
    
    counter <- reactiveVal(1)
    rSelections <- reactiveVal(list())
    rActiveSelectionId <- reactiveVal()
    rActiveSelection <- reactive({
      x <- base::Find(function(elem) { elem$id == rActiveSelectionId() }, rSelections())
      x$rSelection
    }) %>% bindEvent(rActiveSelectionId())
    selectionServer <- sampleTableServer("activeSelectionTable", metadata, rActiveSelection)
    
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
    
    observe({
      updateSelectizeInput(
        session,
        "activeSelectionSelect",
        choices = c(lapply(rSelections(), function(elem) { elem$id }), list("New")),
        selected = rActiveSelectionId()
      )
    }) %>% bindEvent(rSelections())
    
    observe({
      rActiveSelectionId(input$activeSelectionSelect)
    }) %>% bindEvent(input$activeSelectionSelect)
    
    output$sampleSelectionsUi <- renderUI({
      if(is.null(rActiveSelection())) {
        newSelectionUi(ns)
      } else {
        selectionUi(ns)
      }
    })
    
    return(
      list(
        selections = rSelections,
        activeSelectionId = rActiveSelectionId,
        activeSelection = rActiveSelection
      )
    )
  })
}
