#' Create RiboCrypt app
#' @return RiboCrypt shiny app
#' @export

RiboCrypt_app <- function() {
  ui <- fluidPage(
    sidebarLayout(
      sidebarPanel(
        tabsetPanel(id = "tabset",
                    tabPanel("Browser",
                             selectizeInput(
                               inputId = "dff",
                               label = "Select an experiment",
                               choices = "",
                               selected = "",
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
  
  server <- function(input, output, ...) {
    # Initialize experiment list
    experimentList <- reactive(list.experiments()$name)
    # Initialize experiment selector
    observeEvent(experimentList, {
      updateSelectizeInput(
        inputId = "dff",
        choices = experimentList(),
        selected = experimentList()[1]
      )
    })
    # Initialize variables
    df <- reactive(read.experiment(input$dff))
    cds <- reactive({
      dep <- df()
      if (!is.null(dep)) {
        loadRegion(dep) 
      } else { NULL }
    })
    libs <- reactive({
      dep <- df()
      if (!is.null(dep)) {
        bamVarName(dep) 
      } else { NULL }
    })
    # Initialize remaining selectors
    observeEvent(cds, {
      updateSelectizeInput(inputId = 'gene',
                           choices = names(cds()),
                           selected = names(cds())[1],
                           server = TRUE
                           )
    }, ignoreNULL = TRUE)
    observeEvent(libs, {
      updateSelectizeInput(
        inputId = "library",
        choices = libs(),
        selected = libs()[min(length(libs()), 9)]
        )
    }, ignoreNULL = TRUE)
  
    # Main plot
    display_region <- eventReactive(input$go, {
      if (input$gene %in% c("", "NULL")) {
        names(cds()[1])
      } else { input$gene }
    })
    dff <- eventReactive(input$go, {
      libs_to_pick <- if (is.null(input$library)) {
        libs()[1]
      } else { input$library }
      df()[which(libs() %in% libs_to_pick),]
    })
    customRegions <- eventReactive(input$go, {
      if(isTRUE(input$useCustomRegions)) {
        orfs_flt <- fread("~/custom_regions.csv")
        orfs_flt_grl <- GRanges(orfs_flt) %>% groupGRangesBy(.,.$names)
      } else { NULL }
    })
    output$c <- renderPlotly({
      isolate({
        gc()
        stopifnot(is(cds(), "GRangesList"))
        stopifnot(is(libs(), "character"))
        stopifnot(is(df(), "experiment"))
        stopifnot(is(input$gene, "character"))
        })
      RiboCrypt::multiOmicsPlot_ORFikExp(display_range = display_region(),
                                           df = dff(),
                                           display_sequence = "nt",
                                           reads = force(outputLibs(dff(), type = "pshifted", output.mode = "envirlist", naming = "full")),
                                           trailer_extension = input$extendTrailers,
                                           leader_extension = input$extendLeaders,
                                           annotation = "cds",
                                           viewMode = ifelse(input$viewMode, "genomic","tx"),
                                           kmers = input$kmer,
                                           frames_type = input$frames_type,
                                           custom_regions = customRegions())
        })
    
    # Setup for structure viewer
    dynamicVisible <- reactiveVal(FALSE)
    observeEvent(input$selectedRegion, {
      if (!is.null(input$selectedRegion) && !is.null(customRegions())) {
        dynamicVisible(TRUE)
      } else {
        dynamicVisible(FALSE)
      }
    })
    observeEvent(input$dynamicClose, {
      dynamicVisible(FALSE)
    })
    output$dynamic <- renderNGLVieweR({
      paste("~", "sequences", input$selectedRegion, "ranked_0.pdb", sep = "/") %>%
        NGLVieweR() %>%
        stageParameters(backgroundColor = "white") %>% 
        addRepresentation("cartoon")
    })
    output$variableUi <- renderUI({
      if (dynamicVisible()) {
        fluidRow(
          actionButton("dynamicClose", "Close"),
          NGLVieweROutput("dynamic")
        )
        } else {}
      })
  }
  
  shinyApp(ui, server)
}



