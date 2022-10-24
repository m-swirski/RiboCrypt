#' Create RiboCrypt app
#' @param validate.experiments logical, default TRUE, set to FALSE
#' to allow starting the app with malformed experiments, be careful
#' will crash if you try to load that experiment!
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
    ),
    tags$li(
      tags$a(href = "/metadata", "metadata"),
    )
  )
  envirs <- c(browser_env = new.env(), heatmap_env = new.env())

  brochureApp(
    ## Pages
    # Main menu
    landing_page(nav_links),
    # Genome / Transcriptome browser
    browser_page(nav_links, validate.experiments, envirs[["browser_env"]]),
    # Heatmap browser
    heatmap_page(nav_links, validate.experiments, envirs[["heatmap_env"]]),
    #metadata download and display
    metadata_page(nav_links)
  )
}

