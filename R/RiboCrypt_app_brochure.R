#' Create RiboCrypt app
#' @import shiny
#' @importFrom NGLVieweR renderNGLVieweR
#' @return RiboCrypt shiny app
#' @export
#'


RiboCrypt_app_brochure <- function(validate.experiments = TRUE) {
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
    if(!requireNamespace("brochure")) stop("devtools::install_github('ColinFay/brochure')"),
    ## Pages
    # Main menu
    landing_page(nav_links),
    # Genome / Transcriptome browser
    browser_page(nav_links, validate.experiments),
    # Heatmap browser
    heatmap_page(nav_links, validate.experiments)
  )
}

