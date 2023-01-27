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
                   codon_filter_input_select(ns),
                   codon_score_input_select(ns)
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
      # Organism / study objects
      org_and_study_changed_checker(input, output, session)
      # Gene objects
      valid_genes_subset <- reactive(filterTranscripts(df(), stopOnEmpty = FALSE,
                                                       minFiveUTR = 3)) %>%
        bindEvent(rv_changed(), ignoreNULL = TRUE)
      tx <- reactive({loadRegion(df(), part = "mrna", names.keep = valid_genes_subset())}) %>%
        bindEvent(rv_changed(), ignoreNULL = TRUE)
      cds <- reactive(loadRegion(df(), part = "cds", names.keep = valid_genes_subset())) %>%
        bindEvent(rv_changed(), ignoreNULL = TRUE)
      gene_name_list <- reactive(get_gene_name_categories(df())) %>%
        bindCache(rv$curval) %>%
        bindEvent(rv_changed(), ignoreNULL = TRUE)

      # Update main side panels
      all_is_gene <- TRUE
      study_and_gene_observers(input, output, session)

      # Main plot, this code is only run if 'plot' is pressed
      mainPlotControls <- eventReactive(input$go,
                     click_plot_codon_main_controller(input, tx, cds, libs, df))

      coverage <- reactive({
        message("-- Codon analysis: ")
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
          return(dt)
        } else {
          print("This is not a mRNA / valid mRNA")
          return(NULL)
        }
      }) %>%
        bindCache(mainPlotControls()$normalization,
                  ORFik:::name_decider(mainPlotControls()$dff, naming = "full"),
                  mainPlotControls()$filter_value)


      output$c <- renderPlotly({
        message("-- Plotting codon usage")
        score_column <-
        if (input$codon_score == "percentage") {
          coverage()$relative_to_max_score
        } else if (input$codon_score == "dispersion(NB)") {
          coverage()$dispersion_txNorm
        } else if (input$codon_score == "alpha(DMN)") {
          coverage()$alpha
        } else if (input$codon_score == "sum") {
          coverage()$sum
        }
        plotly::ggplotly(ggplot(coverage(),
                                aes(type, seqs, fill = score_column)) +
                           geom_tile(color = "white") +
                           scale_fill_gradient2(low = "blue", high = "orange",
                                                mid = "white") +
                           theme(axis.text.y = element_text(family = "monospace")) +
                           facet_wrap(coverage()$variable, ncol = 4))
      }) %>%
        bindEvent(coverage(), ignoreNULL = TRUE)
    }
  )
}
