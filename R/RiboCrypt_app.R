#' Create RiboCrypt app
#' @param validate.experiments logical, default TRUE, set to FALSE
#' to allow starting the app with malformed experiments, be careful
#' will crash if you try to load that experiment!
#' @param options list of arguments, default
#'  \code{list("launch.browser" = ifelse(interactive(), TRUE, FALSE))}
#' @import shiny bslib ORFik NGLVieweR
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
      browser_ui("browser", validate.experiments = validate.experiments),
      heatmap_ui("heatmap", validate.experiments = validate.experiments),
      metadata_ui("metadata")))

  server <- function(input, output, session) {
    browser_server("browser", all_experiments = all_exp)
    heatmap_server("heatmap", all_experiments = all_exp)
    metadata_server("metadata")
  }

  shinyApp(ui, server, options = options)
}


RiboCrypt_app_modular <- RiboCrypt_app
