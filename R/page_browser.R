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
                   tx_input_select(ns),
                   library_input_select(ns),
                   frame_type_select(ns),
                   sliderInput(ns("kmer"), "K-mer length", min = 1, max = 20, value = 1)
          ),
          tabPanel("Settings",
                   numericInput(ns("extendLeaders"), "5' extension", 0),
                   numericInput(ns("extendTrailers"), "3' extension", 0),
                   checkboxInput(ns("viewMode"), label = "Genomic View", value = FALSE),
                   checkboxInput(ns("useCustomRegions"), label = "Protein structures", value = FALSE),
                   checkboxInput(ns("other_tx"), label = "Full annotation", value = FALSE),
                   checkboxInput(ns("summary_track"), label = "Summary top track", value = FALSE),
                   frame_type_select(ns, "summary_track_type", "Select summary display type"),
                   export_format_of_plot(ns)
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
      # Static values
      genomes <- unique(all_exp$organism)
      experiments <- all_exp$name

      # Set reactive values
      org <- reactive(input$genome)
      rv <- reactiveValues(lstval="",curval="") # Store current and last genome
      rv_changed <- reactiveVal(NULL) # Did genome change?
      df <- reactive(get_exp(input$dff, experiments, env))
      observeEvent(df(), update_rv(rv, df), priority = 2)
      observe(update_rv_changed(rv, rv_changed), priority = 1) %>%
        bindEvent(rv$curval)
      tx <- reactive({loadRegion(df())}) %>%
        bindEvent(rv_changed(), ignoreNULL = T)
      cds <- reactive(loadRegion(df(), part = "cds")) %>%
        bindEvent(rv_changed(), ignoreNULL = T)
      gene_name_list <- reactive(get_gene_name_categories(df())) %>%
        bindEvent(rv_changed(), ignoreNULL = T)
      libs <- reactive(bamVarName(df()))
      # Update main side panels
      observeEvent(org(), experiment_update_select(org, all_exp, experiments))
      observeEvent(gene_name_list(), gene_update_select(gene_name_list))
      observeEvent(input$gene, tx_update_select(isolate(input$gene),
                      gene_name_list), ignoreNULL = TRUE, ignoreInit = T)
      observeEvent(libs(), library_update_select(libs))

      # Main plot controller, this code is only run if 'plot' is pressed
      mainPlotControls <- eventReactive(input$go,
        click_plot_browser_main_controller(input, tx, cds, libs, df))
      # Main plot, this code is only run if 'plot' is pressed
      output$c <- renderPlotly(click_plot_browser(mainPlotControls, session))
      ### NGLVieweR (protein structures) ###
      # TODO: Move as much as possible of protein stuff out of page_browser

      # Setup reactive values needed for structure viewer
      dynamicVisible <- reactiveVal(FALSE)
      selectedRegion <- reactiveVal(NULL)
      selectedRegionProfile <- reactive({
        req(selectedRegion())
        result <- cds()[names(cds()) == selectedRegion()] %>%
          getRiboProfile(mainPlotControls()$reads[[1]]) %>%
          (function (x) { x$count[seq.int(1, length(x$count), 3)] })()
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
      pdb_file_exists <- reactive(pdb_exists(pdb_file))
      output$dynamic <- renderNGLVieweR(
        protein_struct_render(pdb_file_exists, selectedRegionProfile, pdb_file))
      # Variable UI logic
      output$variableUi <- renderUI(
        protein_struct_plot(selectedRegionProfile, dynamicVisible,
                            pdb_file_exists, session))
    }
  )
}
