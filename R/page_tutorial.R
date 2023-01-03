tutorial_ui <- function(id, label = "tutorial") {
  ns <- NS(id)
  tabPanel(title = "tutorial", icon = icon("question"),
    uiOutput(ns('tutorial'))
  )
}


tutorial_server <- function(id) {
  moduleServer(
    id,
    function(input, output, session) {
      output$tutorial <- renderUI({
        tutorial <- HTML(markdown::mark_html(
          knitr::knit(
            system.file("rmd", "tutorial.rmd", package = "RiboCrypt"),
            quiet = TRUE, output = paste0(tempfile(), ".html")),
          options = list(base64_images = TRUE, standalone = FALSE, toc = FALSE)))

        fluidPage(
          fluidRow(
            column(10, tutorial, offset = 1)
          )
        )
      })
    }
  )
}
