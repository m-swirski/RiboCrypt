

landing_page <- function(nav_links) {
  page(
    href = "/",
    ui = function(request) {
      fluidPage(
        h1("Welcome to RiboCrypt!"),
        nav_links
      )
    },
    server = function(input, output, session) {
      # observe({
      #   print(reactiveValuesToList(session$clientData)$url_pathname)
      # })
    }
  )
}

browser_page <- function(nav_links, validate.experiments = TRUE,
                         envir = .GlobalEnv) {
  page(
    href = nav_links$rel_paths["browser"],
    ui = function(request) {
      fluidPage(
        nav_links,
        sidebarLayout(
          sidebarPanel(
            tabsetPanel(id = "tabset",
                        tabPanel("Browser",
                                 experiment_input_select(list.experiments(validate = validate.experiments)$name),
                                 gene_input_select(),
                                 library_input_select(),
                                 frame_type_select(),
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
      df <- reactive(read.experiment(input$dff, output.env = envir))
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
        cds_display <- {
          if (input$gene %in% c("", "NULL")) {
            cds()[0]
          } else {
            if (input$gene %in% names(cds())) {
              cds()[input$gene]
            } else cds()[0]
          }
        }
        reactiveValues(dff = dff,
                       display_region = display_region,
                       customRegions = customRegions,
                       extendTrailers = input$extendTrailers,
                       extendLeaders = input$extendLeaders,
                       viewMode = input$viewMode,
                       kmerLength = input$kmer,
                       frames_type = input$frames_type,
                       cds_display = cds_display)
      })
      output$c <- renderPlotly({
        read_type <- ifelse(dir.exists(file.path(dirname(mainPlotControls()$dff$filepath[1]), "cov_RLE")),
                            "cov", "pshifted")
        message("Using type: ", read_type)
        time_before <- Sys.time()
        a <- RiboCrypt::multiOmicsPlot_ORFikExp(display_range = mainPlotControls()$display_region,
                                           df = mainPlotControls()$dff,
                                           display_sequence = "nt",
                                           reads = force(outputLibs(mainPlotControls()$dff, type = read_type, output.mode = "envirlist", naming = "full", BPPARAM = BiocParallel::SerialParam())),
                                           trailer_extension = mainPlotControls()$extendTrailers,
                                           leader_extension = mainPlotControls()$extendLeaders,
                                           annotation = mainPlotControls()$cds_display, # TODO; THIS IS WRONG; GIVES ERROR; FIX!
                                           viewMode = ifelse(mainPlotControls()$viewMode, "genomic","tx"),
                                           kmers = mainPlotControls()$kmerLength,
                                           frames_type = mainPlotControls()$frames_type,
                                           custom_regions = mainPlotControls()$customRegions)
        cat("lib loading + Coverage calc: "); print(round(Sys.time() - time_before, 2))
        return(a)
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

heatmap_page <- function(nav_links, validate.experiments = TRUE,
                         envir = .GlobalEnv) {
  page(
    href = nav_links$rel_paths["heatmap"],
    ui = function(request) {
      fluidPage(
        nav_links,
        sidebarLayout(
          sidebarPanel(
            tabsetPanel(id = "tabset",
                        tabPanel("heatmap",
                                 experiment_input_select(list.experiments(validate = validate.experiments)$name),
                                 gene_input_select(),
                                 library_input_select(),
                                 selectizeInput(
                                   inputId = "region",
                                   label = "View region",
                                   choices = c("Start codon", "Stop codon"),
                                   selected = "Start codon",
                                   multiple = FALSE
                                 ),
                                 selectizeInput(
                                   inputId = "normalization",
                                   label = "Normalization",
                                   choices = c("transcriptNormalized", "zscore", "sum", "log10sum"),
                                   selected = "transcriptNormalized",
                                   multiple = FALSE
                                 ),
                                 numericInput("readlength_min", "Min Readlength", 26),
                                 numericInput("readlength_max", "Max Readlength", 34)

                        ),
                        tabPanel("Navigate",
                                 numericInput("extendLeaders", "5' extension", 30),
                                 numericInput("extendTrailers", "3' extension", 30)
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
      df <- reactive(read.experiment(input$dff, output.env = envir))
      # TODO: make sure to update valid genes, when 5' and 3' extension is updated!
      valid_genes_subset <- reactive(filterTranscripts(df()))
      tx <- reactive(loadRegion(df(), part = "mrna", names.keep = valid_genes_subset()))
      cds <- reactive(loadRegion(df(), part = "cds", names.keep = valid_genes_subset()))
      libs <- reactive(bamVarName(df()))

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
        display_region <- {
          if (input$gene %in% c("", "NULL")) {
            "all"
          } else { input$gene }
        }
        dff <- {
          libs_to_pick <- if (is.null(input$library)) {
            libs()[1]
          } else { input$library }
          df()[which(libs() %in% libs_to_pick),]
        }

        cds_display <- {
          if (input$gene %in% c("", "NULL")) {
            cds()[0]
          } else {
            if ("all" %in% input$gene) {
              cds()
            } else if (input$gene %in% names(cds())) {
              cds()[input$gene]
            } else cds()[0]
          }
        }
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
        read_type <- ifelse(dir.exists(file.path(dirname(mainPlotControls()$dff$filepath[1]), "cov_RLE_List")),
                            "covl", "pshifted")
        message("-- Using type: ", ifelse(read_type == "pshifted", "ofst", "covl"))
        message("Environment: ")
        print(envExp(mainPlotControls()$dff))
        if (length(mainPlotControls()$cds_display) > 0) {
          print("This is a mRNA")
          time_before <- Sys.time()
          force(outputLibs(mainPlotControls()$dff, type = read_type, output.mode = "envir", naming = "full", BPPARAM = BiocParallel::SerialParam()))
          cat("Library loading: "); print(round(Sys.time() - time_before, 2))
          message("-- Data loading complete")
          # Pick start or stop region
          if (mainPlotControls()$region == "Start codon") {
            region <- groupGRangesBy(startSites(mainPlotControls()$cds_display, TRUE, TRUE, TRUE))
          } else {
            region <- groupGRangesBy(stopSites(mainPlotControls()$cds_display, TRUE, TRUE, TRUE))
          }
          class(region)
          print(class(get(bamVarName(mainPlotControls()$dff, skip.condition = FALSE, skip.replicate = FALSE)[1], envir = envExp(mainPlotControls()$dff))))
          time_before <- Sys.time()
          dt <- windowPerReadLength(region, tx()[names(region)],
                                    reads = get(bamVarName(mainPlotControls()$dff, skip.condition = FALSE, skip.replicate = FALSE)[1], envir = envExp(mainPlotControls()$dff)),
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
          print("This is not a mRNA")
          return(NULL)
        }
      })
    }
  )
}

# dt <- coveragePerTiling(startRegion(cds["ENST00000674920"], mrna["ENST00000674920"], is.sorted = TRUE,
#                                     upstream = 5,
#                                     downstream = 5), reads = get(bamVarName(df, skip.condition = FALSE, skip.replicate = FALSE)[1], envir = envExp(df)), fraction = 28, as.data.table = TRUE)
# dt <- coverageScorings(dt, scoring = "transcriptNormalized")
# ORFik:::coverageHeatMap(dt)

metadata_page <- function(nav_links) {
  page(
    href = nav_links$rel_paths["metadata"],
    ui <- function(request) {
      fluidPage(
        nav_links,
        titlePanel("Metadata search page"),
        sidebarLayout(
          sidebarPanel(
            textInput("accession", "Study accession number (SRP/GEO/PRJNA)")
            ,
            actionButton("go", "Plot"),
          ),
          mainPanel(
            textOutput("abstract"),
            dataTableOutput("metadata")
          )
        )
      )
    },
    server <- function(input, output, session) {
      md <- eventReactive(input$go,{
        abstract <- capture.output({
          meta_dt <- download.SRA.metadata(input$accession)
        })
        reactiveValues(abstract = abstract,
                       meta_dt = meta_dt)
      })
      output$abstract <- renderText(md()$abstract)
      output$metadata <- renderDataTable(md()$meta_dt)

    }

  )
}

