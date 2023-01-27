tutorial_ui <- function() {
  tabPanel(title = "tutorial", icon = icon("question"),
           fluidPage(
             fluidRow(
               column(10, includeHTML(system.file("rmd", "tutorial.html", package = "RiboCrypt")), offset = 1)
             )
           )
  )
}
