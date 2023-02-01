heatmap_ui <- function(id, label = "Heatmap", all_exp) {
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
                   tx_input_select(ns),
                   library_input_select(ns, FALSE),
                   heatmap_region_select(ns),
                   normalization_input_select(ns),
                   numericInput(ns("readlength_min"), "Min Readlength", 26),
                   numericInput(ns("readlength_max"), "Max Readlength", 34)),
          tabPanel("Settings",
                   numericInput(ns("extendLeaders"), "5' extension", 30),
                   numericInput(ns("extendTrailers"), "3' extension", 30),
                   checkboxInput(ns("summary_track"), label = "Summary top track", value = FALSE),
                   checkboxInput(ns("p_shifted"), label = "p_shifted", value = TRUE)), ),
        actionButton(ns("go"), "Plot", icon = icon("rocket")), ),
      mainPanel(
        tabsetPanel(type = "tabs",
                    tabPanel("Heatmap", plotlyOutput(outputId = ns("c"), height = "500px") %>% shinycssloaders::withSpinner(color="#0dc5c1"),
                             uiOutput(ns("variableUi"))),
                    tabPanel("Shift table", tableOutput(outputId = ns("shift_table")) %>% shinycssloaders::withSpinner(color="#0dc5c1"))))
  ))
}

heatmap_server <- function(id, all_experiments, env) {
  moduleServer(
    id,
    function(input, output, session, all_exp = all_experiments) {
      # Organism / study objects
      org_and_study_changed_checker(input, output, session)
      # Gene objects
      valid_genes_subset <- reactive(filterTranscripts(df(), stopOnEmpty = FALSE,
                                          minFiveUTR = 0, minThreeUTR = 0)) %>%
        bindEvent(rv_changed(), ignoreNULL = TRUE)
      length_table <- reactive(optimizedTranscriptLengths(df(), TRUE, TRUE, FALSE)) %>%
        bindEvent(rv_changed(), ignoreNULL = TRUE)
      tx <- reactive({loadRegion(df(), part = "mrna", names.keep = valid_genes_subset())}) %>%
        bindEvent(rv_changed(), ignoreNULL = TRUE)
      cds <- reactive(loadRegion(df(), part = "cds", names.keep = valid_genes_subset())) %>%
        bindEvent(rv_changed(), ignoreNULL = TRUE)
      gene_name_list <- reactive(get_gene_name_categories(df())) %>%
        bindEvent(rv_changed(), ignoreNULL = TRUE)
      libs <- reactive(bamVarName(df()))

      # Update main side panels
      all_is_gene <- TRUE
      study_and_gene_observers(input, output, session)

      # Main plot, this code is only run if 'plot' is pressed
      mainPlotControls <- eventReactive(input$go,
              click_plot_heatmap_main_controller(input, tx, cds, libs, df))

      coverage <- reactive(heatmap_data(mainPlotControls, tx, length_table)) %>%
        bindCache(mainPlotControls()$extendLeaders, mainPlotControls()$extendTrailers,
                  mainPlotControls()$normalization, mainPlotControls()$region,
                  ORFik:::name_decider(mainPlotControls()$dff, naming = "full"),
                  mainPlotControls()$readlength_min,
                  mainPlotControls()$readlength_max)

    output$c <- renderPlotly({
      message("-- Plotting heatmap")
      main_plot <- coverageHeatMap(coverage(), scoring = mainPlotControls()$normalization,
                                   legendPos = "bottom")
      plot_list <- if (mainPlotControls()$summary_track) {
        heights <- c(0.2,0.8)
        list(pSitePlot(coverage(), forHeatmap = TRUE), main_plot)
      } else {
        heights <- 1
        list(main_plot)
      }
      return(subplot(plot_list, nrows = length(plot_list), heights = heights, shareX = TRUE, titleY = TRUE))
    }) %>%
      bindEvent(coverage(), ignoreNULL = TRUE)
    output$shift_table <- renderTable(mainPlotControls()$shift_table)
  }
  )
}
