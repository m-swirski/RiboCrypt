dynamicTabsetPanel <- function(id, selected) {
  function (...) {
    tabsetPanel(..., id = id, selected = selected)
  }
}

sampleSelectionsUi <- function(id) {
  ns <- NS(id)
  uiOutput(ns("sampleSelectionsUi"))
}

nonEmptySelectionUi <- function(id, label) {
  tabPanel(
    label,
    sampleTableUi(id),
    value = id
  )
}

sampleSelectionsServer <- function(id, metadata, rInitialSelection) {
  moduleServer(id, function(input, output, session){
    ns <- session$ns
    
    # list(id = ns("1"), observers = c(), selectedSamples = reactiveVal())
    activeModules <- reactiveVal(list())
    
    counter <- reactiveVal(1)
    
    emptySelectionId <- ns("0")
    emptySelectionUi <- tabPanel(
        "New selection",
        fluidRow(
          column(3, actionButton(ns("saveAsSelection"), "Save selection"))
        ),
        value = emptySelectionId
      )
    
    observe({
      selection <- rInitialSelection()
      modules <- activeModules()
      
      count <- counter()
      counter(count + 1)
      
      newModule <- sampleTableServer(as.character(count), metadata, selection)
      
      activeModules(
        c(modules, list(newModule))
        )
    }) %>% bindEvent(input$saveAsSelection)
    
    output$sampleSelectionsUi <- renderUI({
      emptySelection <- list(emptySelectionUi)
      
      if (length(activeModules()) == 0) {
        tabsetArgs <- emptySelection
        selectedTabId <- emptySelectionId
      } else {
        nonEmptySelections <- lapply(activeModules(), function(x) {
          nonEmptySelectionUi(x$id, paste0("Samples selection ", x$label))
        })
        tabsetArgs <- c(nonEmptySelections, emptySelection)
        selectedTabId <- tail(activeModules(), 1)[[1]]$id
      }
      
      base::do.call(dynamicTabsetPanel(id = "sampleSelectionsTabset", selected = selectedTabId), tabsetArgs)
    })
    
    return(
      list(
        selections = activeModules,
        activeSelectionId = reactive(
          input$sampleSelectionsTabset
        )
      )
    )
  })
}
