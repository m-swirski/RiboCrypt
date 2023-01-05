browser_ui = function(id, label = "Browser", all_exp) {
  ns <- NS(id)
  genomes <- unique(all_exp$organism)
  experiments <- all_exp$name
  tabPanel(
    title = "browser", icon = icon("chart-line"),
    sidebarLayout(
      sidebarPanel(
        tabsetPanel(
          tabPanel("Browser",
                   organism_input_select(c("ALL", genomes), ns),
                   experiment_input_select(experiments, ns),
                   gene_input_select(ns),
                   library_input_select(ns),
                   frame_type_select(ns),
                   sliderInput(ns("kmer"), "K-mer length", min = 1, max = 20, value = 1)
          ),
          tabPanel("Settings",
                   numericInput(ns("extendLeaders"), "5' extension", 0),
                   numericInput(ns("extendTrailers"), "3' extension", 0),
                   checkboxInput(ns("viewMode"), label = "Genomic View", value = FALSE),
                   checkboxInput(ns("useCustomRegions"), label = "Protein structures", value = FALSE),
                   checkboxInput(ns("summary_track"), label = "Summary top track", value = FALSE),
                   frame_type_select(ns, "summary_track_type", "Select summary display type")
          ),
        ),
        actionButton(ns("go"), "Plot", icon = icon("rocket")),
      ),
      mainPanel(
        plotlyOutput(outputId = ns("c")) %>% shinycssloaders::withSpinner(color="#0dc5c1"),
        uiOutput(ns("variableUi"))
      )
    )
  )
}

browser_server <- function(id, all_experiments, env) {
  moduleServer(
    id,
    function(input, output, session, all_exp = all_experiments) {
      # Loading selected experiment and related data
      genomes <- unique(all_exp$organism)
      experiments <- all_exp$name
      # Set reactive values
      org <- reactive(input$genome)
      df <- reactive(read.experiment(input$dff, output.env = env)) #, output.env = envir))
      tx <- reactive(loadRegion(df()))
      cds <- reactive(loadRegion(df(), part = "cds"))
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

      # Main plot, this code is only run if 'plot' is pressed
      mainPlotControls <- eventReactive(input$go, {
        display_region <- observed_gene(isolate(input$gene), tx)
        cds_display <- observed_cds(isolate(input$gene),cds)
        dff <- observed_exp_subset(isolate(input$library), libs, df)
        customRegions <- load_custom_regions(isolate(input$useCustomRegions), df)

        reads <- load_reads(dff, "cov")
        reactiveValues(dff = dff,
                       display_region = display_region,
                       customRegions = customRegions,
                       extendTrailers = input$extendTrailers,
                       extendLeaders = input$extendLeaders,
                       summary_track = input$summary_track,
                       summary_track_type = input$summary_track_type,
                       viewMode = input$viewMode,
                       kmerLength = input$kmer,
                       frames_type = input$frames_type,
                       cds_display = cds_display,
                       reads = reads)
      })

      output$c <- renderPlotly({
        time_before <- Sys.time()
        a <- RiboCrypt::multiOmicsPlot_ORFikExp(
          display_range = mainPlotControls()$display_region,
          df = mainPlotControls()$dff,
          display_sequence = "nt",
          reads = mainPlotControls()$reads,
          trailer_extension = mainPlotControls()$extendTrailers,
          leader_extension = mainPlotControls()$extendLeaders,
          #annotation = mainPlotControls()$cds_display,
          viewMode = ifelse(mainPlotControls()$viewMode, "genomic","tx"),
          kmers = mainPlotControls()$kmerLength,
          frames_type = mainPlotControls()$frames_type,
          custom_regions = mainPlotControls()$customRegions,
          input_id = session$ns("selectedRegion"),
          summary_track = mainPlotControls()$summary_track,
          summary_track_type = mainPlotControls()$summary_track_type
          )
        cat("lib loading + Coverage calc: "); print(round(Sys.time() - time_before, 2))
        return(a)
      })

      ### NGLVieweR ###
      # TODO: Move as much as possible of protein stuff out of page_browser
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
