browser_ui = function(id,  all_exp, browser_options, gene_names_init,
                      libs, label = "Browser") {
  ns <- NS(id)
  genomes <- unique(all_exp$organism)
  experiments <- all_exp$name
  init_tx <- gene_names_init[label == browser_options["default_gene"],]
  tabPanel(
    title = "browser", icon = icon("chart-line"),
    sidebarLayout(
      jqui_resizable(sidebarPanel(
        tabsetPanel(
          tabPanel("Browser",
                   organism_input_select(c("ALL", genomes), ns),
                   experiment_input_select(experiments, ns, browser_options),
                   gene_input_select(ns, FALSE, browser_options),
                   tx_input_select(ns, FALSE, init_tx),
                   library_input_select(ns, TRUE, libs),
                   frame_type_select(ns),
                   sliderInput(ns("kmer"), "K-mer length", min = 1, max = 20, value = 1),
                   helper_button_redirect_call()
          ),
          tabPanel("Settings",
                   numericInput(ns("extendLeaders"), "5' extension", 0),
                   numericInput(ns("extendTrailers"), "3' extension", 0),
                   checkboxInput(ns("viewMode"), label = "Genomic View", value = FALSE),
                   checkboxInput(ns("useCustomRegions"), label = "Protein structures", value = FALSE),
                   checkboxInput(ns("other_tx"), label = "Full annotation", value = FALSE),
                   checkboxInput(ns("add_uorfs"), label = "Add uORFs", value = FALSE),
                   checkboxInput(ns("summary_track"), label = "Summary top track", value = FALSE),
                   frame_type_select(ns, "summary_track_type", "Select summary display type"),
                   export_format_of_plot(ns)
          )
        ),
         actionButton(ns("go"), "Plot", icon = icon("rocket")), width=3
      )),
      mainPanel(
        jqui_resizable(plotlyOutput(outputId = ns("c"), height = "500px")) %>% shinycssloaders::withSpinner(color="#0dc5c1"),
        uiOutput(ns("variableUi")
      ),width=9)
    )
  )
}

browser_server <- function(id, all_experiments, env, df, experiments,
                           tx, cds, libs, org, gene_name_list, rv,
                           browser_options) {
  moduleServer(
    id,
    function(input, output, session, all_exp = all_experiments) {
      # browser()
      study_and_gene_observers(input, output, session)
      # Main plot controller, this code is only run if 'plot' is pressed
      mainPlotControls <- eventReactive(input$go,
        click_plot_browser_main_controller(input, tx, cds, libs, df),
        ignoreInit = !as.logical(browser_options["plot_on_start"]),
        ignoreNULL = FALSE)
      # Main plot, this code is only run if 'plot' is pressed
      output$c <- renderPlotly(click_plot_browser(mainPlotControls, session)) %>%
        bindEvent(mainPlotControls(),
                  ignoreInit = FALSE,
                  ignoreNULL = TRUE)
      # Protein display
      module_protein(input, output, session)
      return(rv)
    }
  )
}
