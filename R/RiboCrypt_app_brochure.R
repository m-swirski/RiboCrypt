#' Create RiboCrypt app
#' @param validate.experiments logical, default TRUE, set to FALSE
#' to allow starting the app with malformed experiments, be careful
#' will crash if you try to load that experiment!
#' @import shiny
#' @importFrom NGLVieweR renderNGLVieweR
#' @return RiboCrypt shiny app
#' @export
#' @examples
#'
#' ## To run in RSTUDIO server using ssh
#' ## A proxy url path is made, so we need to assign that
#' ## First run the app, and lcopy url part after port
#' ## should look something like this: /p/3b3a7b68/
#'
RiboCrypt_app_brochure <- function(validate.experiments = TRUE,
                                   options = list("launch.browser" = ifelse(interactive(), TRUE, FALSE)),
                                   base.url = ".") {
  if(!requireNamespace("brochure")) {
    message("To run the brochure app, please install brochure using: ")
    stop("devtools::install_github('ColinFay/brochure')")
  } else library(brochure)

  nav_links <- nav_links_creator(base.url)
  envirs <- c(browser_env = new.env(), heatmap_env = new.env())

  brochureApp(
    ## Pages
    # INIT page
    start_page(nav_links),
    # Main menu
    landing_page(nav_links),
    # Genome / Transcriptome browser
    browser_page(nav_links, validate.experiments, envirs[["browser_env"]]),
    # Heatmap browser
    heatmap_page(nav_links, validate.experiments, envirs[["heatmap_env"]]),
    #metadata download and display
    metadata_page(nav_links),
    options = options
  )
}




nav_links_creator <- function(b = "/") {
  if (!is(b, "character")) stop("base.url must be character!")
  library(shiny)
  paths <- c(home = b,
             browser = paste0(b, "browser"),
             heatmap = paste0(b, "heatmap"),
             metadata = paste0(b, "metadata"))
  nav_links <- tags$ul(
    tags$li(tags$a(href = paths["home"], "home"),),
    tags$li(tags$a(href = paths["browser"], "browser"), ),
    tags$li(tags$a(href = paths["heatmap"], "heatmap"), ),
    tags$li(tags$a(href = paths["metadata"], "metadata"),))
  nav_links$rel_paths <- paths

  return(nav_links)
}

RiboCrypt_app_brochure_old <- function(validate.experiments = TRUE) {
  if(!requireNamespace("brochure")) {
    message("To run the brochure app, please install brochure using: ")
    stop("devtools::install_github('ColinFay/brochure')")
  } else library(brochure)

  library(shiny)
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
    RiboCrypt:::landing_page(nav_links),
    # Genome / Transcriptome browser
    RiboCrypt:::browser_page(nav_links, validate.experiments, envirs[["browser_env"]]),
    # Heatmap browser
    RiboCrypt:::heatmap_page(nav_links, validate.experiments, envirs[["heatmap_env"]]),
    #metadata download and display
    RiboCrypt:::metadata_page(nav_links)
  )
}

