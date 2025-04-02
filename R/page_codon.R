codon_ui <- function(id, all_exp, browser_options, libs, label = "Codon") {
  ns <- NS(id)
  genomes <- unique(all_exp$organism)
  experiments <- all_exp$name
  tabPanel(
    title = "Codon dwell time", icon = icon("layer-group"),
    sidebarLayout(
      sidebarPanel(
        tabsetPanel(
          tabPanel("Codon",
                   organism_input_select(c("ALL", genomes), ns),
                   experiment_input_select(experiments, ns, browser_options),
                   gene_input_select(ns),
                   tx_input_select(ns),
                   library_input_select(ns, TRUE, libs),
                   codon_filter_input_select(ns),
                   codon_score_input_select(ns),
                   checkboxInput(ns("differential"), label = "Differential", value = FALSE),
                   )),
        plot_button(ns("go"))
      ),
      mainPanel(
        plotlyOutput(outputId = ns("c"), height = "750px") %>%
          shinycssloaders::withSpinner(color="#0dc5c1")))
  )
}

codon_server <- function(id, all_experiments, env, df, experiments, tx, cds,
                         libs, org, gene_name_list, rv) {
  moduleServer(
    id,
    function(input, output, session, all_exp = all_experiments) {
      length_table <- reactive(optimizedTranscriptLengths(df(), TRUE, TRUE)) %>%
        bindCache(rv$curval) %>%
        bindEvent(rv$changed, ignoreNULL = TRUE)
      # Update main side panels
      all_is_gene <- TRUE
      study_and_gene_observers(input, output, session)

      # Main plot, this code is only run if 'plot' is pressed
      mainPlotControls <- eventReactive(input$go,
                     click_plot_codon_main_controller(input, tx, cds, libs, df,
                                                      length_table))

      coverage <- reactive({
        message("-- Codon analysis: ")
        if (length(mainPlotControls()$cds_display) > 0) {
          print("Valid input")
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
      output$c <- renderPlotly(click_plot_codon(input, coverage)) %>%
        bindEvent(coverage(), ignoreInit = FALSE, ignoreNULL = TRUE)
      return(rv)
    }
  )
}
