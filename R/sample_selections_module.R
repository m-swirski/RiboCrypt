sampleSelectionsUi <- function(id) {
  ns <- NS(id)
  uiOutput(ns("sampleSelectionsUi"))
}


sampleSelectionsServer <- function(id, metadata, rInitialSelectedSamples) {
  moduleServer(id, function(input, output, session){
    ns <- session$ns
    counter <- reactiveVal(1)
    selectedModule <- ""
    # list(id = ns("1"), observers = c(), selectedSamples = reactiveVal())
    activeModules <- reactiveVal(c())
    
    observe({
      currentSelection <- rInitialSelectedSamples()
      currentActiveModules <- activeModules()
      
      count <- counter()
      counter(count + 1)
      
      newModule <- sampleTableServer(as.character(count), metadata, currentSelection)
      
      activeModules(base::append(currentActiveModules, newModule))
    }) %>% bindEvent(input$saveAsSelection)
    
    output$sampleSelectionsUi <- renderUI({
      nonEmptySelections <- lapply(activeModules()$id, function(id) {
        nonEmptySelectionUi(ns(id), paste0("Samples ", id))
      })
      emptySelection <- list(emptySelectionUi(ns))
      
      if (length(nonEmptySelections) == 0) {
        tabsetArgs <- emptySelection
      } else {
        tabsetArgs <- c(nonEmptySelections, emptySelection)
      }
      
      base::do.call(tabsetPanel, tabsetArgs)
    })
    
    return(c())
  })
}

emptySelectionUi <- function(ns) {
  tabPanel(
    "New selection",
      fluidRow(
        column(3, actionButton(ns("saveAsSelection"), "Save selection"))
      )
    )
}

nonEmptySelectionUi <- function(id, label) {
  tabPanel(
    label,
    sampleTableUi(id)
  )
}