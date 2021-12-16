
#' Create RiboCrypt app
#' @return RiboCrypt shiny app
#' @export

RiboCrypt_app <- function() {
  ui <- fluidPage(
    sidebarLayout(
      sidebarPanel(
        tabsetPanel(id = "tabset",
                    tabPanel("browser",
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
                    )
                    ,
                             tabPanel("Navigate",
                             numericInput("extendLeaders",
                               "5' extension",
                               0
                             ), numericInput("extendTrailers",
                                             "3' extension",
                                             0),
                             checkboxInput("viewMode", label = "Genomic View", value = FALSE),
                                             ),

                    ),

        actionButton("go", "Plot")
      ),
      mainPanel(
        plotlyOutput(outputId = "c")
      )
    )
  )

  server <- function(input, output, ...) {
    print("Start")
    # Init variables
    exp.list <- list.experiments()$name
    d <- read.experiment(exp.list[1])
    last.ex <- reactiveVal(exp.list[1])
    cd <- loadRegion(d)
    lib <- bamVarName(d)
    # Init reactive variables
    df <- reactiveVal(d)
    cds <- reactiveVal(cd)
    libs <- reactiveVal(lib)
    print("Init done variables")
    # Init input boxes
    updateSelectizeInput(inputId = 'gene',
                         choices = names(cd),
                         selected = names(cd[1]),
                         server = TRUE
    )
    updateSelectizeInput(
      inputId = "library",
      choices = lib,
      selected = lib[min(length(lib), 9)]
    )
    updateSelectizeInput(
      inputId = "dff",
      choices = exp.list,
      selected = exp.list[1],
    )
    print("Init done Selectize")

    v <- reactiveValues(doPlot = FALSE)

    observeEvent(input$go, {
      # 0 will be coerced to FALSE
      # 1+ will be coerced to TRUE
      v$doPlot <- input$go
    })



    # Plot and updates
    output$c <- renderPlotly({
      if (v$doPlot == FALSE) return()
      isolate ({
        gc()
        stopifnot(is(cds(), "GRangesList"))
        stopifnot(is(libs(), "character"))
        stopifnot(is(df(), "experiment"))
        stopifnot(is(input$gene, "character"))
        print("Inside plotly")
        paste("Input gene", input$gene)
        # Init display and lib at start run, else it crash now
        display_region <- if (input$gene %in% c("", "NULL")) {
          names(cds()[1])
        } else input$gene
        libs_to_pick <- if (is.null(input$library)) {
          libs()[1]
        } else input$library

        dff <- df()[which(libs() %in% libs_to_pick),]
        show(dff)
        observeEvent(input$dff, {
          cat("Changed experiment")
          if (input$dff != "" & input$dff != last.ex()) {
            print(input$dff)
            df(read.experiment(input$dff))
            last.ex(input$dff)
            cds(loadRegion(df()))
            libs(bamVarName(df()))

            updateSelectizeInput(
              inputId = 'gene',
              choices = names(cds()),
              selected = names(cds()[1]),
              server = TRUE
            )
            updateSelectizeInput(
              inputId = "library",
              choices = libs(),
              selected = libs()[min(length(lib), 9)]
            )
          }
        }, ignoreInit = TRUE)

        RiboCrypt::multiOmicsPlot_ORFikExp(display_range = display_region,
                                           df = dff,
                                           display_sequence = "nt",
                                           reads = force(outputLibs(dff, type = "pshifted", output.mode = "envirlist", naming = "full")),
                                           trailer_extension = input$extendTrailers,
                                           leader_extension = input$extendLeaders,
                                           annotation = "cds",
                                           viewMode = ifelse(input$viewMode, "genomic","tx"),
                                           kmers = input$kmer,
                                           frames_type = input$frames_type)
      })
    })
  }


  shinyApp(ui, server)
}


