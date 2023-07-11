browser_allsamp_ui = function(id,  all_exp, browser_options,
                      metadata, label = "Browser_allsamp") {
  ns <- NS(id)
  genomes <- unique(all_exp$organism)
  experiments <- all_exp$name
  tabPanel(
    title = "browser", icon = icon("chart-line"),
    sidebarLayout(
      jqui_resizable(sidebarPanel(
        tabsetPanel(
          tabPanel("Browser",
                   organism_input_select(c("ALL", genomes), ns),
                   experiment_input_select(experiments, ns, browser_options),
                   gene_input_select(ns, FALSE, browser_options),
                   metadata_input_select(ns, metadata = metadata),
                   sliderInput(ns("kmer"), "K-mer length", min = 1, max = 20,
                               value = as.numeric(browser_options["default_kmer"])),
                   helper_button_redirect_call()
          ),
          tabPanel("Settings",
                   numericInput(ns("extendLeaders"), "5' extension", 0),
                   numericInput(ns("extendTrailers"), "3' extension", 0),
                   textInput(ns("customSequence"), label = "Custom sequences highlight", value = NULL),
                   tx_input_select(ns, FALSE),
                   checkboxInput(ns("viewMode"), label = "Genomic View", value = FALSE),
                   checkboxInput(ns("useCustomRegions"), label = "Protein structures", value = FALSE),
                   checkboxInput(ns("other_tx"), label = "Full annotation", value = FALSE),
                   checkboxInput(ns("add_uorfs"), label = "uORF annotation", value = FALSE),
                   checkboxInput(ns("expression_plot"), label = "Add expression plot", value = FALSE),
                   checkboxInput(ns("summary_track"), label = "Summary top track", value = FALSE),
                   frame_type_select(ns, "summary_track_type", "Select summary display type"),
                   export_format_of_plot(ns)
          )
        ),
        actionButton(ns("go"), "Plot", icon = icon("rocket")), width=3
      )),
      mainPanel(
        jqui_resizable(plotOutput(outputId = ns("c"), height = "500px")) %>% shinycssloaders::withSpinner(color="#0dc5c1"), width=9)
    )
  )
}

browser_allsamp_server <- function(id, all_experiments, env, df, experiments,
                           tx, cds, libs, org, gene_name_list, rv, metadata) {
  moduleServer(
    id,
    function(input, output, session, all_exp = all_experiments) {
      # browser()
      study_and_gene_observers(input, output, session)
      # Main plot controller, this code is only run if 'plot' is pressed
      mainPlotControls <- eventReactive(input$go,
                                        click_plot_browser_allsamp_controller(input, tx, libs, df),
                                        ignoreInit = TRUE,
                                        ignoreNULL = FALSE)
      # Main plot, this code is only run if 'plot' is pressed
      output$c <- renderPlot(click_plot_browser_allsamples(mainPlotControls)) %>%
        bindCache(ORFik:::name_decider(mainPlotControls()$dff, naming = "full"),
                  input$tx, input$other_tx, input$add_uorfs,
                  input$extendTrailers, input$extendLeaders,
                  input$plot_export_format, input$genomic_region,
                  input$summary_track, input$summary_track_type,
                  input$viewMode, input$kmer, input$frames_type,input$customSequence) %>%
        bindEvent(mainPlotControls(),
                  ignoreInit = FALSE,
                  ignoreNULL = TRUE)
    }
  )
}
