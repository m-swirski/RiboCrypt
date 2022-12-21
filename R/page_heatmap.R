heatmap_ui <- function(id, label = "Heatmap", validate.experiments = T,
                       all_exp = list.experiments(validate = validate.experiments)) {
  ns <- NS(id)
  genomes <- unique(all_exp$organism)
  experiments <- all_exp$name
  tabPanel(
    title = "heatmap", icon = icon("layer-group"),
    sidebarLayout(
      sidebarPanel(
        tabsetPanel(
          tabPanel("Heatmap",
                   organism_input_select(c("ALL", genomes), ns),
                   experiment_input_select(experiments, ns),
                   gene_input_select(ns),
                   library_input_select(ns),
                   heatmap_region_select(ns),
                   normalization_input_select(ns),
                   numericInput(ns("readlength_min"), "Min Readlength", 26),
                   numericInput(ns("readlength_max"), "Max Readlength", 34)),
          tabPanel("Navigate",
                   numericInput(ns("extendLeaders"), "5' extension", 30),
                   numericInput(ns("extendTrailers"), "3' extension", 30)), ),
        actionButton(ns("go"), "Plot", icon = icon("rocket")), ),
      mainPanel(
        plotlyOutput(outputId = ns("c")),
        uiOutput(ns("variableUi"))))
  )
}

heatmap_server <- function(id, all_experiments) {
  moduleServer(
    id,
    function(input, output, session, all_exp = all_experiments) {
      # Loading selected experiment and related data
      genomes <- unique(all_exp$organism)
      experiments <- all_exp$name
      # Set reactive values
      org <- reactive(input$genome)
      df <- reactive(read.experiment(input$dff)) #, output.env = envir))
      # TODO: make sure to update valid genes, when 5' and 3' extension is updated!
      valid_genes_subset <- reactive(filterTranscripts(df(), stopOnEmpty = FALSE, minThreeUTR = 0))
      tx <- reactive(loadRegion(df(), part = "mrna", names.keep = valid_genes_subset()))
      cds <- reactive(loadRegion(df(), part = "cds", names.keep = valid_genes_subset()))
      libs <- reactive(bamVarName(df()))

      observeEvent(org(), {
        orgs_safe <- if (isolate(org()) == "ALL") {
          unique(all_exp$organism)
        } else isolate(org())
        picks <- experiments[all_exp$organism %in% orgs_safe]
        updateSelectizeInput(
          inputId = "dff",
          choices = picks,
          selected = picks[1],
          server = FALSE
        )
      })

      # Gene selector
      observeEvent(tx(), {
        updateSelectizeInput(
          inputId = 'gene',
          choices = c("all", names(tx())),
          selected = "all",
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

      # Main plot, this code is only run if 'plot' is pressed
      mainPlotControls <- eventReactive(input$go, {
        display_region <- observed_gene_heatmap(isolate(input$gene), tx)
        cds_display <- observed_cds_heatmap(isolate(input$gene),cds)
        dff <- observed_exp_subset(isolate(input$library), libs, df)
        reactiveValues(dff = dff,
                       display_region = display_region,
                       extendTrailers = input$extendTrailers,
                       extendLeaders = input$extendLeaders,
                       viewMode = input$viewMode,
                       cds_display = cds_display,
                       region = input$region,
                       readlength_min = input$readlength_min,
                       readlength_max = input$readlength_max,
                       normalization = input$normalization)
      })
      output$c <- renderPlotly({
        message("-- Plot region: ", mainPlotControls()$region)
        filepath1 <- mainPlotControls()$dff$filepath[1]
        read_type <- ifelse(dir.exists(file.path(dirname(filepath1), "cov_RLE_List")),
                            "covl", "pshifted")
        message("-- Using type: ", ifelse(read_type == "pshifted", "ofst", "covl"))
        message("Environment: ")
        print(envExp(mainPlotControls()$dff))
        if (length(mainPlotControls()$cds_display) > 0) {
          print("This is a mRNA")
          time_before <- Sys.time()
          force(outputLibs(
            mainPlotControls()$dff,
            type = read_type,
            output.mode = "envir",
            naming = "full",
            BPPARAM = BiocParallel::SerialParam()))
          cat("Library loading: "); print(round(Sys.time() - time_before, 2))
          message("-- Data loading complete")
          # Pick start or stop region
          region <- observed_cds_point(mainPlotControls)
          time_before <- Sys.time()
          dt <- windowPerReadLength(region, tx()[names(region)],
                                    reads = get(bamVarName(mainPlotControls()$dff, FALSE, FALSE, FALSE, FALSE)[1], envir = envExp(mainPlotControls()$dff)),
                                    pShifted = FALSE, upstream = mainPlotControls()$extendLeaders,
                                    downstream = mainPlotControls()$extendTrailers - 1,
                                    scoring = mainPlotControls()$normalization,
                                    acceptedLengths = seq(mainPlotControls()$readlength_min, mainPlotControls()$readlength_max),
                                    drop.zero.dt = TRUE, append.zeroes = TRUE)
          print(paste("Number of rows in dt:", nrow(dt)))
          cat("Coverage calc: "); print(round(Sys.time() - time_before, 2))
          return(subplot(list(coverageHeatMap(dt, scoring = mainPlotControls()$normalization,
                                              legendPos = "bottom"))))
        } else {
          print("This is not a mRNA / valid mRNA")
          return(NULL)
        }
      })
    }
  )
}
