codon_ui <- function(id, label = "Codon", all_exp) {
  ns <- NS(id)
  genomes <- unique(all_exp$organism)
  experiments <- all_exp$name
  tabPanel(
    title = "codon", icon = icon("layer-group"),
    sidebarLayout(
      sidebarPanel(
        tabsetPanel(
          tabPanel("Codon",
                   organism_input_select(c("ALL", genomes), ns),
                   experiment_input_select(experiments, ns),
                   gene_input_select(ns),
                   tx_input_select(ns),
                   library_input_select(ns),
                   codon_filter_input_select(ns)
                   )),
        actionButton(ns("go"), "Plot", icon = icon("rocket")), ),
      mainPanel(
        plotlyOutput(outputId = ns("c"), height = "750px") %>%
          shinycssloaders::withSpinner(color="#0dc5c1")))
  )
}

codon_server <- function(id, all_experiments, env) {
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
                                                       minFiveUTR = 3)) %>%
        bindEvent(rv_changed(), ignoreNULL = T)
      tx <- reactive({loadRegion(df(), part = "mrna", names.keep = valid_genes_subset())}) %>%
        bindEvent(rv_changed(), ignoreNULL = T)
      cds <- reactive(loadRegion(df(), part = "cds", names.keep = valid_genes_subset())) %>%
        bindEvent(rv_changed(), ignoreNULL = T)
      gene_name_list <- reactive(get_gene_name_categories(df())) %>%
        bindEvent(rv_changed(), ignoreNULL = T)
      libs <- reactive(bamVarName(df()))
      # Update main side panels
      observeEvent(org(), experiment_update_select(org, all_exp, experiments))
      observeEvent(gene_name_list(), gene_update_select_heatmap(gene_name_list))
      observeEvent(input$gene, tx_update_select(isolate(input$gene),
                                                gene_name_list, "all"), ignoreNULL = TRUE, ignoreInit = TRUE)
      observeEvent(libs(), library_update_select(libs))

      # Main plot, this code is only run if 'plot' is pressed
      mainPlotControls <- eventReactive(input$go,
                     click_plot_codon_main_controller(input, tx, cds, libs, df))

      output$c <- renderPlotly({
        message("-- Plot codons: ")
        if (length(mainPlotControls()$cds_display) > 0) {
          print("This is a mRNA")
          filter_val <- mainPlotControls()$filter_value
          print(paste("Filter value:", filter_val))
          print(class(mainPlotControls()$reads[[1]]))
          time_before <- Sys.time()
          dt <- codon_usage_exp(mainPlotControls()$dff,
                                reads = mainPlotControls()$reads,
                                cds = mainPlotControls()$cds_display,
                                mrna = tx()[names(mainPlotControls()$cds_display)],
                                min_counts_cds_filter = filter_val)
          print(paste("Number of rows in dt:", nrow(dt)))
          cat("Coverage calc: "); print(round(Sys.time() - time_before, 2))
          return(subplot(list(codon_usage_plot(dt,
                              score_column = dt$relative_to_max_score))))
        } else {
          print("This is not a mRNA / valid mRNA")
          return(NULL)
        }
      })
    }
  )
}
