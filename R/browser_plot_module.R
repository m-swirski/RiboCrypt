browserPlotUi <- function(id, gene_names_init, browser_options) {
  ns <- NS(id)
  all_isoforms <- subset(gene_names_init, label == browser_options["default_gene"])
  div(
    fluidRow(
      column(2, gene_input_select(ns, FALSE, browser_options)),
      column(2, tx_input_select(ns, FALSE, all_isoforms)),
      column(1, NULL, plot_button(ns("go")))
    ),
    fluidRow(
      plotlyOutput(ns("plot"))
    )
  )
}

browserPlotServer <- function(id, tx, cds, libs, df, browser_options, rSelectedSamples) {
  
  createCoverageSubplot <- function(coverage) {
    # TODO 
  }
  
  createAnnotationSubplot <- function() {
    # TODO
  }
  
  createSequenceSubplot <- function() {
    # TODO
  }
  
  moduleServer(id, function(input, output, session) {
    introns_width <- as.numeric(browser_options["collapsed_introns_width"])
    # Main plot controller, this code is only run if 'plot' is pressed
    mainPlotControls <- reactive({
      req(rSelectedSamples())
      click_plot_browser_main_controller(
        tx, 
        cds, 
        libs, 
        df,
        selectedGene = input$gene, 
        selectedTx = input$tx, 
        otherTx = FALSE, 
        addUorfs = FALSE,
        addTranslons = FALSE,
        collapsedIntrons = FALSE,
        collapsedIntronsWidth = introns_width,
        genomicRegion = "",
        extendLeaders = 0,
        extendTrailers = 0,
        zoomRange = "",
        viewMode = FALSE,
        selectedLibraries = c(""),
        withFrames = TRUE,
        exportFormat = c("svg", "png")[1],
        summaryTrack = FALSE,
        summaryTrackType = c("lines", "columns", "stacks", "area", "heatmap", "animate")[1],
        kmer = 9,
        framesType = c("lines", "columns", "stacks", "area", "heatmap", "animate")[1],
        framesSubset = "all",
        normalization = normalizations("metabrowser")[1],
        customSequence = "",
        logScale = FALSE,
        phyloP = FALSE,
        mapability = FALSE,
        expressionPlot = FALSE
      )
    }) %>% bindEvent(input$go, ignoreInit = TRUE, ignoreNULL = FALSE)
    
    bottom_panel <- reactive(bottom_panel_shiny(mainPlotControls))  %>%
      bindCache(mainPlotControls()$hash_bottom) %>%
      bindEvent(mainPlotControls(), ignoreInit = FALSE, ignoreNULL = TRUE)
    
    browser_plot <- reactive(browser_track_panel_shiny(mainPlotControls, bottom_panel(), session)) %>%
      bindCache(mainPlotControls()$hash_browser) %>%
      bindEvent(bottom_panel(), ignoreInit = FALSE, ignoreNULL = TRUE)
    
    output$plot <- renderPlotly(browser_plot()) %>%
      bindCache(mainPlotControls()$hash_browser) %>%
      bindEvent(browser_plot(), ignoreInit = FALSE, ignoreNULL = TRUE)
    
  })
}