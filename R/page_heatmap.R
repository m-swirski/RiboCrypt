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
                   numericInput(ns("extendTrailers"), "3' extension", 30)), ),
        actionButton(ns("go"), "Plot", icon = icon("rocket")), ),
      mainPanel(
        plotlyOutput(outputId = ns("c")) %>% shinycssloaders::withSpinner(color="#0dc5c1"),
        uiOutput(ns("variableUi"))))
  )
}

heatmap_server <- function(id, all_experiments, env) {
  moduleServer(
    id,
    function(input, output, session, all_exp = all_experiments) {
      # Loading selected experiment and related data
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
      valid_genes_subset <- reactive(filterTranscripts(df(), stopOnEmpty = FALSE,
                                                       minThreeUTR = 0)) %>%
        bindEvent(rv_changed(), ignoreNULL = TRUE)
      tx <- reactive({loadRegion(df(), part = "mrna", names.keep = valid_genes_subset())}) %>%
        bindEvent(rv_changed(), ignoreNULL = TRUE)
      cds <- reactive(loadRegion(df(), part = "cds", names.keep = valid_genes_subset())) %>%
        bindEvent(rv_changed(), ignoreNULL = TRUE)
      gene_name_list <- reactive(get_gene_name_categories(df())) %>%
        bindEvent(rv_changed(), ignoreNULL = TRUE)
      libs <- reactive(bamVarName(df()))

      # Update main side panels
      observeEvent(org(), experiment_update_select(org, all_exp, experiments))
      observeEvent(gene_name_list(), gene_update_select_heatmap(gene_name_list))
      observeEvent(input$gene, tx_update_select(isolate(input$gene),
              gene_name_list, "all"), ignoreNULL = TRUE, ignoreInit = TRUE)
      observeEvent(libs(), library_update_select(libs))

      # Main plot, this code is only run if 'plot' is pressed
      mainPlotControls <- eventReactive(input$go,
              click_plot_heatmap_main_controller(input, tx, cds, libs, df))

      output$c <- renderPlotly({
        message("-- Plot region: ", mainPlotControls()$region)
        if (length(mainPlotControls()$cds_display) > 0) {
          print("This is a mRNA")
          print(class(mainPlotControls()$reads[[1]]))
          # Pick start or stop region
          region <- observed_cds_point(mainPlotControls)
          time_before <- Sys.time()
          dt <- windowPerReadLength(region, tx()[names(region)],
                                    reads = mainPlotControls()$reads[[1]],
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
