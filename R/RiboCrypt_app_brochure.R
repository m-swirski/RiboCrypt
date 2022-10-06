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
    )
  )

  landing_page <- function() {
    page(
      href = "/",
      ui = function(request) {
        fluidPage(
          h1("Welcome to RiboCrypt!"),
          nav_links
        )
      },
      server = function(input, output, session) {

      }
    )
  }

  browser_page <- function() {
    page(
      href = "/browser",
      ui = function(request) {
        fluidPage(
          nav_links,
          sidebarLayout(
            sidebarPanel(
              tabsetPanel(id = "tabset",
                          tabPanel("Browser",
                                   selectizeInput(
                                     inputId = "dff",
                                     label = "Select an experiment",
                                     choices = list.experiments()$name,
                                     selected = head(list.experiments()$name),
                                     multiple = FALSE
                                   ),
                                   selectizeInput(
                                     inputId = "gene",
                                     choices = "",
                                     selected = "",
                                     label = "Select a gene",
                                     multiple = FALSE
                                   ),
                                   selectizeInput(
                                     inputId = "library",
                                     label = "Select libraries",
                                     choices = "",
                                     selected = "",
                                     multiple = TRUE
                                   ),
                                   selectizeInput(
                                     inputId = "frames_type",
                                     label = "Select frames display type",
                                     choices = c("lines", "columns", "stacks", "area"),
                                     selected = "lines",
                                     multiple = FALSE
                                   ),
                                   sliderInput("kmer", "K-mer length", min = 1, max = 20, value = 1)
                          ),
                          tabPanel("Navigate",
                                   numericInput("extendLeaders", "5' extension", 0),
                                   numericInput("extendTrailers", "3' extension", 0),
                                   checkboxInput("viewMode", label = "Genomic View", value = FALSE),
                                   checkboxInput("useCustomRegions", label = "Use custom regions", value = FALSE)
                          ),
              ),
              actionButton("go", "Plot"),
            ),
            mainPanel(
              plotlyOutput(outputId = "c"),
              uiOutput("variableUi")
            )
          )
        )
      },
      server <- function(input, output, ...) {
        # Loading selected experiment and related data
        df <- reactive(read.experiment(input$dff))
        tx <- reactive(loadRegion(df()))
        cds <- reactive(loadRegion(df(), part = "cds"))
        libs <- reactive(bamVarName(df()))

        # Gene selector
        observeEvent(tx(), {
          updateSelectizeInput(
            inputId = 'gene',
            choices = names(tx()),
            selected = names(tx())[1],
            server = TRUE
          )
        })

        # Library selector
        observeEvent(libs(), {
          updateSelectizeInput(
            inputId = "library",
            choices = libs(),
            selected = libs()[min(length(libs()), 9)]
          )
        })

        # Main plot
        mainPlotControls <- eventReactive(input$go, {
          display_region <- {
            if (input$gene %in% c("", "NULL")) {
              names(tx()[1])
            } else { input$gene }
          }
          dff <- {
            libs_to_pick <- if (is.null(input$library)) {
              libs()[1]
            } else { input$library }
            df()[which(libs() %in% libs_to_pick),]
          }
          customRegions <- {
            if(isTRUE(input$useCustomRegions)) {
              orfs_flt <- fread("~/custom_regions.csv")
              orfs_flt_grl <- GRanges(orfs_flt) %>% groupGRangesBy(.,.$names)
            } else { NULL }
          }
          reactiveValues(dff = dff,
                         display_region = display_region,
                         customRegions = customRegions,
                         extendTrailers = input$extendTrailers,
                         extendLeaders = input$extendLeaders,
                         viewMode = input$viewMode,
                         kmerLength = input$kmer,
                         frames_type = input$frames_type)
        })
        output$c <- renderPlotly({
          read_type <- ifelse(dir.exists(file.path(dirname(mainPlotControls()$dff$filepath[1]), "cov_RLE")),
                              "cov", "pshifted")
          message("Using type: ", read_type)
          RiboCrypt::multiOmicsPlot_ORFikExp(display_range = mainPlotControls()$display_region,
                                             df = mainPlotControls()$dff,
                                             display_sequence = "nt",
                                             reads = force(outputLibs(mainPlotControls()$dff, type = read_type, output.mode = "envirlist", naming = "full", BPPARAM = BiocParallel::SerialParam())),
                                             trailer_extension = mainPlotControls()$extendTrailers,
                                             leader_extension = mainPlotControls()$extendLeaders,
                                             annotation = cds(), # TODO; THIS IS WRONG; GIVES ERROR; FIX!
                                             viewMode = ifelse(mainPlotControls()$viewMode, "genomic","tx"),
                                             kmers = mainPlotControls()$kmerLength,
                                             frames_type = mainPlotControls()$frames_type,
                                             custom_regions = mainPlotControls()$customRegions)
        })

        # Setup variables needed for structure viewer
        selectedRegion <- reactiveVal(NULL)
        dynamicVisible <- reactiveVal(FALSE)
        # When user clicks on region
        # start displaying structure viewer
        # and set selected structure to one which was clicked
        observeEvent(input$selectedRegion, {
          if (!is.null(input$selectedRegion) && !is.null(mainPlotControls()$customRegions)) {
            selectedRegion(input$selectedRegion)
            dynamicVisible(TRUE)
          } else {
            selectedRegion(NULL)
            dynamicVisible(FALSE)
          }
        })
        # When user clicks close button
        # stop displaying structure viewer
        # and set selected structure to NULL
        observeEvent(input$dynamicClose, {
          selectedRegion(NULL)
          dynamicVisible(FALSE)
        })
        # NGL viewer widget
        output$dynamic <- renderNGLVieweR({
          paste("~", "sequences", selectedRegion(), "ranked_0.pdb", sep = "/") %>%
            NGLVieweR() %>%
            stageParameters(backgroundColor = "white") %>%
            addRepresentation("cartoon")
        })
        # Variable UI logic
        output$variableUi <- renderUI({
          if (dynamicVisible()) {
            fluidRow(
              actionButton("dynamicClose", "Close"),
              NGLVieweROutput("dynamic")
            )
          } else {}
        })
      }
    )
  }

  brochureApp(
    if(!requireNamespace("brochure")) stop("devtools::install_github('ColinFey/brochure'"),
    # Pages
    landing_page(),
    browser_page()
  )

}

