#' Create RiboCrypt app
#' @param validate.experiments logical, default TRUE, set to FALSE
#' to allow starting the app with malformed experiments, be careful
#' will crash if you try to load that experiment!
#' @param options list of arguments, default
#'  \code{list("launch.browser" = ifelse(interactive(), TRUE, FALSE))}
#' @import shiny bslib ORFik NGLVieweR
#' @importFrom shinycssloaders withSpinner
#' @importFrom markdown mark_html
#' @importFrom knitr knit
#' @return RiboCrypt shiny app
#' @export
#' @examples
#'
#' ## To run in RSTUDIO server using ssh
#' ## A proxy url path is made, so we need to assign that
#' ## First run the app, and lcopy url part after port
#' ## should look something like this: /p/3b3a7b68/
#'
RiboCrypt_app <- function(
    validate.experiments = TRUE,
    options = list("launch.browser" = ifelse(interactive(), TRUE, FALSE))) {
  # Load dataset status
  all_exp <- list.experiments(validate = validate.experiments)
  browser_env <- new.env()
  heatmap_env <- new.env()
  addResourcePath(prefix = "images",
                  directoryPath = system.file("images", package = "RiboCrypt"))

  ui <- tagList(
    tags$head(
      tags$link(rel = "icon",
                href = file.path("images", "favicon.ico"),
                type = "image/x-icon")),
    navbarPage(
      lang = "en",
      windowTitle = "RiboCrypt",
      title = withTags(
        a(img(src = file.path("images", "logo.png"),
              alt = "RiboCrypt",
              height = 60))),
      theme = bslib::bs_theme(
        version = 5,
        primary = "#6dbaff", secondary = "#ff7e7e",
        success = "#c0ffa4", font_scale = 1.2, bootswatch = "zephyr"),
      browser_ui("browser", all_exp = all_exp),
      heatmap_ui("heatmap", all_exp = all_exp),
      metadata_ui("metadata"),
      tabPanel(
        title = "tutorial", icon = icon("question"),
        uiOutput('tutorial'))
      )
    )

  server <- function(input, output, session) {
    browser_server("browser", all_exp, browser_env)
    heatmap_server("heatmap", all_exp, heatmap_env)
    metadata_server("metadata")

    output$tutorial <- renderUI({
      tutorial <- HTML(markdown::mark_html(
        knitr::knit(
          system.file("rmd", "tutorial.rmd", package = "RiboCrypt"),
          quiet = TRUE),
        options = list(base64_images = T, standalone = F, toc = F)))

      fluidPage(
        fluidRow(
          column(10, tutorial, offset = 1)
        )
      )
    })
  }

  shinyApp(ui, server, options = options)
}


RiboCrypt_app_modular <- RiboCrypt_app
