quality_ui <- function(id, label = "Periodicity", all_exp) {
  ns <- NS(id)
  genomes <- unique(all_exp$organism)
  experiments <- all_exp$name
  tabPanel(
    title = "periodicity", icon = icon("layer-group"),
    sidebarLayout(
      sidebarPanel(
        tabsetPanel(
          tabPanel("Periodicity",
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

quality_server <- function(id, all_experiments, env) {
  moduleServer(
    id,
    function(input, output, session, all_exp = all_experiments) {
      # Organism / study objects
      org_and_study_changed_checker(input, output, session)
      # Gene objects
      valid_genes_subset <- reactive(filterTranscripts(df(), stopOnEmpty = FALSE,
                                                       minThreeUTR = 0)) %>%
        bindEvent(rv_changed(), ignoreNULL = TRUE)
      tx <- reactive({loadRegion(df(), part = "mrna", names.keep = valid_genes_subset())}) %>%
        bindEvent(rv_changed(), ignoreNULL = TRUE)
      cds <- reactive(loadRegion(df(), part = "cds", names.keep = valid_genes_subset())) %>%
        bindEvent(rv_changed(), ignoreNULL = TRUE)
      gene_name_list <- reactive(get_gene_name_categories(df())) %>%
        bindEvent(rv_changed(), ignoreNULL = TRUE)

      # Update main side panels
      all_is_gene <- TRUE
      study_and_gene_observers(input, output, session)

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
          dt[, frame := as.factor(position %% 3), by = fraction]

          return(periodicity_plot(dt, fft = TRUE))
        } else {
          print("This is not a mRNA / valid mRNA")
          return(NULL)
        }
      })
    }
  )
}
