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

browserPlotServer <- function(id, browser_options, rExperiment, rSelectedSamples) {
  shiny::moduleServer(id, function(input, output, session) {
    introns_width <- as.numeric(browser_options["collapsed_introns_width"])
    loadedExperiment <- shiny::reactive({
      shiny::req(rExperiment())
      shiny::req(rExperiment() != "")
      ORFik::read.experiment(rExperiment(), validate = FALSE)
    })
    tx <- shiny::reactive({
      shiny::req(loadedExperiment())
      ORFik::loadRegion(loadedExperiment())
    }) %>%
      shiny::bindEvent(loadedExperiment())
    cds <- shiny::reactive({
      shiny::req(loadedExperiment())
      ORFik::loadRegion(loadedExperiment(), "cds")
    }) %>%
      shiny::bindEvent(loadedExperiment())
    libs <- shiny::reactive({
      shiny::req(loadedExperiment())
      ORFik::bamVarName(loadedExperiment())
    }) %>%
      shiny::bindEvent(loadedExperiment())
    gene_name_list <- shiny::reactive({
      shiny::req(loadedExperiment())
      get_gene_name_categories(loadedExperiment())
    })
    samples <- shiny::reactive({
      result <- lapply(rSelectedSamples()$filteredSelections, function(x) {
        unlist(x$sample)
      })
      names(result) <- NULL
      result
    })
    shiny::observe({
      shiny::req(gene_name_list())
      gene_update_select(gene_name_list)
    }) %>% shiny::bindEvent(gene_name_list())
    shiny::observe({
      shiny::req(loadedExperiment())
      shiny::req(input$gene)
      shiny::req(input$gene != "")
      tx_update_select(gene = input$gene, gene_name_list = gene_name_list)
    }) %>% shiny::bindEvent(input$gene)
    # Main plot controller, this code is only run if 'plot' is pressed
    mainPlotControls <- shiny::reactive({
      shiny::req(rSelectedSamples())
      click_plot_browser_main_controller(
        tx,
        cds,
        libs,
        loadedExperiment,
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
        framesType = c("area", "lines", "columns", "stacks", "heatmap", "animate")[1],
        framesSubset = NULL,
        normalization = normalizations("metabrowser")[1],
        customSequence = "",
        logScale = FALSE,
        phyloP = FALSE,
        mapability = FALSE,
        expressionPlot = FALSE,
        useFST = TRUE,
        selectedSamples = samples()
      )
    }) %>% shiny::bindEvent(input$go, ignoreInit = TRUE, ignoreNULL = FALSE)

    bottom_panel <- shiny::reactive(bottom_panel_shiny(mainPlotControls)) %>%
      shiny::bindEvent(mainPlotControls(), ignoreInit = FALSE, ignoreNULL = TRUE)

    browser_plot <- shiny::reactive(browser_track_panel_shiny(mainPlotControls, bottom_panel(), session)) %>%
      shiny::bindEvent(bottom_panel(), ignoreInit = FALSE, ignoreNULL = TRUE)

    output$plot <- renderPlotly(browser_plot()) %>%
      shiny::bindEvent(browser_plot(), ignoreInit = FALSE, ignoreNULL = TRUE)
  })
}
