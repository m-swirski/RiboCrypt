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
#' @param browser_options named character vector of browser specific arguments
#' @param init_tab_focus character, default "browser". Which tab to open on
#' init.
#' @import shiny bslib ORFik NGLVieweR ggplot2
#' @importFrom shinycssloaders withSpinner
#' @importFrom markdown mark_html
#' @importFrom shinyjqui jqui_resizable jqui_draggable
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
    all_exp = list.experiments(validate = validate.experiments),
    browser_options = c(), init_tab_focus = "browser") {
  time_before <- Sys.time()

  stopifnot(is(all_exp, "data.table"))
  stopifnot(!is.null(all_exp$name))
  stopifnot(nrow(all_exp) > 0)
  # Set environments
  with_readlengths_env <- new.env()
  without_readlengths_env <- new.env()
  #with_cigar_env <- new.env() # Not used for now
  # Add resource directories
  addResourcePath(prefix = "images",
                  directoryPath = system.file("images", package = "RiboCrypt"))
  addResourcePath(prefix = "rmd",
                  directoryPath = system.file("rmd", package = "RiboCrypt"))
  # Setup variables
  if (!isTruthy(browser_options["default_experiment"])) {
    browser_options["default_experiment"] <- all_exp$name[1]
  }
  if (is.na(browser_options["plot_on_start"])) {
    browser_options["plot_on_start"] <- FALSE
  }
  exp_init <- read.experiment(browser_options["default_experiment"],
                              validate = FALSE)
  names_init <- get_gene_name_categories(exp_init)
  if (!isTruthy(browser_options["default_gene"])) {
    if (!isTruthy(browser_options["default_gene"])) {
      browser_options["default_gene"] <- names_init$label[1]
    }
    stopifnot(browser_options["default_gene"] %in% names_init$label)
  }
  if (!isTruthy(browser_options["default_kmer"])) {
    browser_options["default_kmer"] <- 1
  } else {
    stopifnot(!is.na(as.numeric(browser_options["default_kmer"])))
  }
  if (!isTruthy(browser_options["default_frame_type"])) {
    browser_options["default_frame_type"] <- "lines"
  } else {
    stopifnot(is.character(browser_options["default_frame_type"]))
  }
  libs <- bamVarName(exp_init)

  # User interface
  ui <- tagList(
    rc_header_image(),
    navbarPage(
      lang = "en",
      windowTitle = "RiboCrypt",
      title = rc_title(),
      theme = rc_theme(),
      selected = init_tab_focus,
      browser_ui("browser", all_exp, browser_options, names_init, libs),
      analysis_ui("analysis", all_exp, browser_options, libs),
      metadata_ui("metadata", all_exp),
      tutorial_ui()
    )
  )

  server <- function(input, output, session) {
    org_and_study_changed_checker(input, output, session)

    rv <- browser_server("browser", all_exp, without_readlengths_env, df,
                         experiments, tx, cds, libs, org, gene_name_list, rv,
                         browser_options)
    rv <- analysis_server("analysis", all_exp, without_readlengths_env,
            with_readlengths_env, df, df_with, experiments, tx, cds, libs, org,
            gene_name_list, rv)
    metadata_server("metadata", all_exp)
    cat("Server: "); print(round(Sys.time() - time_before, 2))
  }
  cat("Init: "); print(round(Sys.time() - time_before, 2))
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

