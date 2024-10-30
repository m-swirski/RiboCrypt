tutorial_ui <- function(id) {
  ns <- NS(id)
  tabPanel(title = "tutorial", icon = icon("question"),
           fluidPage(
             htmlOutput(ns("tutorial"))
           )
  )
}

tutorial_server <- function(id) {
  moduleServer(id, function(input, output, session) {
      # browser()
      addResourcePath("rmd", system.file("rmd",package = "RiboCrypt"))
      output$tutorial <- renderUI({
        tags$iframe(
          seamless="seamless",
          src="rmd/tutorial.html",
          style='width:100vw;height:100vh;')
      })
    })
}
