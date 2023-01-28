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
      observeEvent(org(), experiment_update_select(org, all_exp, experiments))
      observeEvent(gene_name_list(), gene_update_select_heatmap(gene_name_list))
      observeEvent(input$gene, tx_update_select(isolate(input$gene),
              gene_name_list, "all"), ignoreNULL = TRUE, ignoreInit = TRUE)
      observeEvent(libs(), library_update_select(libs))

      # Main plot, this code is only run if 'plot' is pressed
      mainPlotControls <- eventReactive(input$go,
              click_plot_heatmap_main_controller(input, tx, cds, libs, df))

      coverage <- reactive({
        message("-- Region: ", mainPlotControls()$region)
        if (length(mainPlotControls()$cds_display) > 0) {
          print("This is a mRNA")
          print(class(mainPlotControls()$reads[[1]]))
          # Pick start or stop region
          point <- observed_cds_point(mainPlotControls)
          windows <- startRegion(point, tx()[names(point)], TRUE,
                                 upstream = mainPlotControls()$extendLeaders,
                                 downstream = mainPlotControls()$extendTrailers - 1)
          length_table_sub <- length_table()[tx_name %in% names(point),]
          if (mainPlotControls()$region == "Start codon") {
            windows <- extend_needed(windows, length_table_sub$utr5_len,
                                     mainPlotControls()$extendLeaders, "up")
            windows <- extend_needed(windows, length_table_sub$cds_len,
                                     mainPlotControls()$extendTrailers - 1, "down")
          } else {
            windows <- extend_needed(windows, length_table_sub$cds_len,
                                     mainPlotControls()$extendLeaders, "up")
            windows <- extend_needed(windows, length_table_sub$utr3_len,
                                     mainPlotControls()$extendTrailers  - 1, "down")
          }

          time_before <- Sys.time()
          # browser()

          dt <- windowPerReadLength(point, tx(),
                                    reads = mainPlotControls()$reads[[1]],
                                    pShifted = FALSE, upstream = mainPlotControls()$extendLeaders,
                                    downstream = mainPlotControls()$extendTrailers - 1,
                                    scoring = mainPlotControls()$normalization,
                                    acceptedLengths = seq(mainPlotControls()$readlength_min, mainPlotControls()$readlength_max),
                                    drop.zero.dt = TRUE, append.zeroes = TRUE,
                                    windows = windows)
          if (!mainPlotControls()$p_shifted){
            
            sdt <- mainPlotControls()$shift_table
            colnames(sdt)[1] <- "readlength"
            dt[, position := position + sdt[readlength == fraction]$offsets_start, by = fraction]
            dt <- dt[position %between% c(- mainPlotControls()$extendLeaders + max(abs(sdt$offsets_start)), 
                                          mainPlotControls()$extendTrailers - 1 - max(abs(sdt$offsets_start)))]
          } 
          print(paste("Number of rows in dt:", nrow(dt)))
          cat("Coverage calc: "); print(round(Sys.time() - time_before, 2))
          return(dt)
        } else {
          print("This is not a mRNA / valid mRNA")
          return(NULL)
        }
      }) %>%
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
