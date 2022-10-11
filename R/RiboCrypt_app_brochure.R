#' Create RiboCrypt app
#' @import shiny
#' @importFrom NGLVieweR renderNGLVieweR
#' @return RiboCrypt shiny app
#' @export
#'


RiboCrypt_app_brochure <- function(validate.experiments = TRUE) {
  if(!requireNamespace("brochure")) {
    message("To run the brochure app, please install brochure using: ")
    stop("devtools::install_github('ColinFay/brochure')")
  } else library(brochure)

  nav_links <- tags$ul(
    tags$li(
      tags$a(href = "/", "home"),
    ),
    tags$li(
      tags$a(href = "/browser", "browser"),
    ),
    tags$li(
      tags$a(href = "/heatmap", "heatmap"),
    )
  )
  brochureApp(
    ## Pages
    # Main menu
    landing_page(nav_links),
    # Genome / Transcriptome browser
    browser_page(nav_links, validate.experiments),
    # Heatmap browser
    heatmap_page(nav_links, validate.experiments)
  )
}

