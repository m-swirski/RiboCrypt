#' Create RiboCrypt app
#' @import shiny
#' @importFrom NGLVieweR renderNGLVieweR
#' @return RiboCrypt shiny app
#' @export
#'


RiboCrypt_app_brochure <- function() {
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
    if(!requireNamespace("brochure")) stop("devtools::install_github('ColinFey/brochure'"),
    ## Pages
    # Main menu
    landing_page(nav_links),
    # Genome / Transcriptome browser
    browser_page(nav_links),
    # Heatmap browser
    heatmap_page(nav_links)
  )
}

