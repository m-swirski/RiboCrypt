#' Create RiboCrypt app
#' @param validate.experiments logical, default TRUE, set to FALSE
#' to allow starting the app with malformed experiments, be careful
#' will crash if you try to load that experiment!
#' @param options list of arguments, default
#'  \code{list("launch.browser" = ifelse(interactive(), TRUE, FALSE))}
#' @param all_exp a data.table, default:
#' \code{list.experiments(validate = validate.experiments)}. Which experiments
#' do you want to allow your app to see, default is all in your system config
#' path.
#' @import shiny bslib ORFik NGLVieweR ggplot2
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
    options = list("launch.browser" = ifelse(interactive(), TRUE, FALSE)),
    all_exp = list.experiments(validate = validate.experiments)) {
  # Set environments
  with_readlengths_env <- new.env()
  without_readlengths_env <- new.env()
  #with_cigar_env <- new.env() # Not used for now
  addResourcePath(prefix = "images",
                  directoryPath = system.file("images", package = "RiboCrypt"))
  addResourcePath(prefix = "rmd",
                  directoryPath = system.file("rmd", package = "RiboCrypt"))

  ui <- tagList(
    rc_header_image(),
    navbarPage(
      lang = "en",
      windowTitle = "RiboCrypt",
      title = rc_title(),
      theme = rc_theme(),
      browser_ui("browser", all_exp = all_exp),
      # heatmap_ui("heatmap", all_exp = all_exp),
      analysis_ui("analysis", all_exp = all_exp),
      metadata_ui("metadata"),
      tutorial_ui()
    )
  )

  server <- function(input, output, session) {
    browser_server("browser", all_exp, without_readlengths_env)
    # heatmap_server("heatmap", all_exp, with_readlengths_env)
    analysis_server("analysis", all_exp, without_readlengths_env,
                    with_readlengths_env)
    metadata_server("metadata")
    # tutorial_server("tutorial")
  }

  shinyApp(ui, server, options = options)
}

RiboCrypt_app_modular <- RiboCrypt_app


rc_header_image <- function() {
  tags$head(
    tags$link(rel = "icon",
              href = file.path("images", "favicon.ico"),
              type = "image/x-icon"))
}

rc_theme <- function() {
  bslib::bs_theme(
    version = 5,
    primary = "#6dbaff", secondary = "#ff7e7e",
    success = "#c0ffa4", font_scale = 1.2, bootswatch = "zephyr")
}

rc_title <- function() {
  withTags(
    a(img(src = file.path("images", "logo.png"),
          alt = "RiboCrypt",
          height = 60)))
}

