browser_ui = function(id, label = "Browser", validate.experiments = T) {
  ns <- NS(id)
  experiments <- list.experiments(validate = validate.experiments)$name
  tabPanel(
    title = "browser", icon = icon("chart-line"),
    sidebarLayout(
      sidebarPanel(
        tabsetPanel(
          tabPanel("Browser",
                   selectizeInput(
                     inputId = ns("dff"),
                     label = "Select an experiment",
                     choices = experiments,
                     selected = head(experiments),
                     multiple = FALSE),
                   selectizeInput(
                     inputId = ns("gene"),
                     choices = "",
                     selected = "",
                     label = "Select a gene",
                     multiple = FALSE),
                   selectizeInput(
                     inputId = ns("library"),
                     label = "Select libraries",
                     choices = "",
                     selected = "",
                     multiple = TRUE),
                   selectizeInput(
                     inputId = ns("frames_type"),
                     label = "Select frames display type",
                     choices = c("lines", "columns", "stacks", "area"),
                     selected = "lines",
                     multiple = FALSE),
                   sliderInput(ns("kmer"), "K-mer length", min = 1, max = 20, value = 1)
          ),
          tabPanel("Navigate",
                   numericInput(ns("extendLeaders"), "5' extension", 0),
                   numericInput(ns("extendTrailers"), "3' extension", 0),
                   checkboxInput(ns("viewMode"), label = "Genomic View", value = FALSE),
                   checkboxInput(ns("useCustomRegions"), label = "Use custom regions", value = FALSE)
          ),
        ),
        actionButton(ns("go"), "Plot", icon = icon("rocket")),
      ),
      mainPanel(
        plotlyOutput(outputId = ns("c")),
        uiOutput(ns("variableUi"))
      )
    )
  )
}

browser_server <- function(id) {
  moduleServer(
    id,
    function(input, output, session) {
      # Loading selected experiment and related data
      df <- reactive(read.experiment(input$dff)) #, output.env = envir))
      tx <- reactive(loadRegion(df()))
      cds <- reactive(loadRegion(df(), part = "cds"))
      libs <- reactive(bamVarName(df()))

      # Gene selector
      observeEvent(tx(), {
        updateSelectizeInput(
          inputId = "gene",
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
            protein_structure_path <- file.path(dirname(df()@fafile), "protein_structure_predictions", "custom_regions.csv")
            if (file.exists(protein_structure_path)) {
              orfs_flt <- fread(protein_structure_path)
              orfs_flt_grl <- GRanges(orfs_flt) %>% groupGRangesBy(.,.$names)
            } else NULL
          } else NULL
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
        read_type <- {
          ifelse(dir.exists(file.path(dirname(dff$filepath[1]), "cov_RLE")), "cov", "pshifted")
        }
        reads <- {
          force(
            outputLibs(
              dff,
              type = read_type,
              output.mode = "envirlist",
              naming = "full",
              BPPARAM = BiocParallel::SerialParam()
            )
          )
        }
        reactiveValues(dff = dff,
                       display_region = display_region,
                       customRegions = customRegions,
                       extendTrailers = input$extendTrailers,
                       extendLeaders = input$extendLeaders,
                       viewMode = input$viewMode,
                       kmerLength = input$kmer,
                       frames_type = input$frames_type,
                       cds_display = cds_display,
                       reads = reads)
      })

      output$c <- renderPlotly({
        filepath1 <- mainPlotControls()$dff$filepath[1]
        read_type <- ifelse(dir.exists(
          file.path(dirname(filepath1), "cov_RLE")),
          "cov", "pshifted")
        message("Using type: ", read_type)

        time_before <- Sys.time()
        a <- RiboCrypt::multiOmicsPlot_ORFikExp(
          display_range = mainPlotControls()$display_region,
          df = mainPlotControls()$dff,
          display_sequence = "nt",
          reads = mainPlotControls()$reads,
          trailer_extension = mainPlotControls()$extendTrailers,
          leader_extension = mainPlotControls()$extendLeaders,
          annotation = mainPlotControls()$cds_display, # TODO; THIS IS WRONG; GIVES ERROR; FIX!
          viewMode = ifelse(mainPlotControls()$viewMode, "genomic","tx"),
          kmers = mainPlotControls()$kmerLength,
          frames_type = mainPlotControls()$frames_type,
          custom_regions = mainPlotControls()$customRegions,
          input_id = session$ns("selectedRegion"))
        cat("lib loading + Coverage calc: "); print(round(Sys.time() - time_before, 2))
        return(a)
      })

      ### NGLVieweR ###

      # Utilities or computing coloring scheme based on Ribo-Seq results
      valuesToColors <- function(vals) {
        palette <- colorRampPalette(c("blue", "green", "yellow", "orange", "red"), bias = 0.5)(1001)
        palette[(vals / max(vals) * 1000) + 1]
      }

      # Setup reactive values needed for structure viewer
      dynamicVisible <- reactiveVal(FALSE)
      selectedRegion <- reactiveVal(NULL)
      selectedRegionProfile <- reactive({
        req(selectedRegion())
        result <- cds()[names(cds()) == selectedRegion()] %>%
          getRiboProfile(mainPlotControls()$reads[[1]]) %>%
          (function (x) { x$count[seq(1, length(x$count), 3)] })()
      })

      # When user clicks on region
      # start displaying structure viewer
      # and set selected structure to one which was clicked
      observeEvent(input$selectedRegion, {
        req(input$selectedRegion)
        selectedRegion(input$selectedRegion)
        dynamicVisible(TRUE)
      })
      # When user clicks close button
      # stop displaying structure viewer
      # and set selected structure to NULL
      observeEvent(input$dynamicClose, {
        selectedRegion(NULL)
        dynamicVisible(FALSE)
      })
      # NGL viewer widget
      protein_structure_dir <- reactive({
        file.path(dirname(df()@fafile), "protein_structure_predictions")
      })
      region_dir <- reactive({
        file.path(protein_structure_dir(), selectedRegion())
      })
      pdb_file <- reactive({
        file.path(region_dir(), "ranked_0.pdb")
      })
      pdb_file_exists <- reactive({
        req(pdb_file())
        file.exists(pdb_file())
      })
      output$dynamic <- renderNGLVieweR({
        req(pdb_file_exists(), selectedRegionProfile())
        pdb_file() %>% NGLVieweR() %>%
          stageParameters(backgroundColor = "white") %>%
          onRender(fetchJS("sequence_viewer_coloring.js"), valuesToColors(selectedRegionProfile()))
      })
      # Variable UI logic
      output$variableUi <- renderUI({
        ns <- session$ns
        req(dynamicVisible(), pdb_file_exists(), selectedRegionProfile())
        fluidRow(
          actionButton(ns("dynamicClose"), "Close"),
          NGLVieweROutput(ns("dynamic"))
        )
      })
    }
  )
}

metadata_ui <- function(id, label = "metadata") {
  ns <- NS(id)
  tabPanel(title = "metadata", icon = icon("rectangle-list"),
           h2("Metadata search page"),
           sidebarLayout(
             sidebarPanel(
               textInput(ns("accession"), "Study accession number (SRP/GEO/PRJNA)"),
               actionButton(ns("go"), "Plot", icon = icon("rocket")),
             ),
             mainPanel(
               textOutput(ns("abstract")),
               dataTableOutput(ns("metadata"))
             )
           )
  )
}

metadata_server <- function(id) {
  moduleServer(
    id,
    function(input, output, session) {
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

heatmap_ui <- function(id, label = "Heatmap", validate.experiments = T) {
  ns <- NS(id)
  experiments <- list.experiments(validate = validate.experiments)$name
  tabPanel(
    title = "heatmap", icon = icon("layer-group"),
    sidebarLayout(
      sidebarPanel(
        tabsetPanel(
          tabPanel("Heatmap",
                   selectizeInput(
                     inputId = ns("dff"),
                     label = "Select an experiment",
                     choices = experiments,
                     selected = head(experiments),
                     multiple = FALSE),
                   selectizeInput(
                     inputId = ns("gene"),
                     choices = "",
                     selected = "",
                     label = "Select a gene",
                     multiple = FALSE),
                   selectizeInput(
                     inputId = ns("library"),
                     label = "Select libraries",
                     choices = "",
                     selected = "",
                     multiple = TRUE),
                   selectizeInput(
                     inputId = ns("region"),
                     label = "View region",
                     choices = c("Start codon", "Stop codon"),
                     selected = "Start codon",
                     multiple = FALSE
                   ),
                   selectizeInput(
                     inputId = ns("normalization"),
                     label = "Normalization",
                     choices =
                       c("transcriptNormalized", "zscore", "sum", "log10sum"),
                     selected = "transcriptNormalized",
                     multiple = FALSE
                   ),
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

heatmap_server <- function(id) {
  moduleServer(
    id,
    function(input, output, session) {
      # Loading selected experiment and related data
      df <- reactive(read.experiment(input$dff)) #, output.env = envir))
      # TODO: make sure to update valid genes, when 5' and 3' extension is updated!
      valid_genes_subset <- reactive(filterTranscripts(df(), stopOnEmpty = FALSE, minThreeUTR = 0))
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
          if (mainPlotControls()$region == "Start codon") {
            region <- groupGRangesBy(
              startSites(mainPlotControls()$cds_display, TRUE, TRUE, TRUE))
          } else {
            region <- groupGRangesBy(
              stopSites(mainPlotControls()$cds_display, TRUE, TRUE, TRUE))
          }
          class(region)
          #print(class(get(bamVarName(mainPlotControls()$dff, skip.condition = FALSE, skip.replicate = FALSE)[1], envir = envExp(mainPlotControls()$dff))))
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

#' Create RiboCrypt app
#' @param validate.experiments logical, default TRUE, set to FALSE
#' to allow starting the app with malformed experiments, be careful
#' will crash if you try to load that experiment!
#' @param options list of arguments, default
#'  \code{list("launch.browser" = ifelse(interactive(), TRUE, FALSE))}
#' @import shiny bslib ORFik NGLVieweR
#' @return RiboCrypt shiny app
#' @export
#' @examples
#'
#' ## To run in RSTUDIO server using ssh
#' ## A proxy url path is made, so we need to assign that
#' ## First run the app, and lcopy url part after port
#' ## should look something like this: /p/3b3a7b68/
#'
RiboCrypt_app_modular <- function(
    validate.experiments = TRUE,
    options = list("launch.browser" = ifelse(interactive(), TRUE, FALSE))) {

  ui <- navbarPage(
    lang = "en",
    title = "RiboCrypt",
    theme = bslib::bs_theme(
      version = 5,
      primary = "#6dbaff", secondary = "#ff7e7e",
      success = "#c0ffa4", font_scale = 1.2, bootswatch = "zephyr"),
    browser_ui("browser", validate.experiments = validate.experiments),
    heatmap_ui("heatmap", validate.experiments = validate.experiments),
    metadata_ui("metadata")
  )

  server <- function(input, output, session) {
    browser_server("browser")
    heatmap_server("heatmap")
    metadata_server("metadata")
  }

  shinyApp(ui, server, options = options)
}
