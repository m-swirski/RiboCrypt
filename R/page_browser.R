browser_ui = function(id,  all_exp, browser_options, gene_names_init,
                      libs, label = "Browser") {
  ns <- NS(id)
  genomes <- unique(all_exp$organism)
  experiments <- all_exp$name
  init_tx <- gene_names_init[label == browser_options["default_gene"],]
  init_libs <- unlist(strsplit(browser_options["default_libs"], "\\|"))
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
                   library_input_select(ns, TRUE, libs, init_libs),
                   frame_type_select(ns, selected =
                                       browser_options["default_frame_type"]),
                   sliderInput(ns("kmer"), "K-mer length", min = 1, max = 20,
                          value = as.numeric(browser_options["default_kmer"])),
                   helper_button_redirect_call(), # Here are settings ->
                   shinyjs::useShinyjs(),
                   rclipboardSetup(),
                   tags$head(
                     tags$style(HTML('#clip{background-color:orange}'))
                   )
          ),
          tabPanel("Settings",
                   numericInput(ns("extendLeaders"), "5' extension", 0),
                   numericInput(ns("extendTrailers"), "3' extension", 0),
                   textInput(ns("customSequence"), label = "Custom sequences highlight", value = NULL),
                   textInput(ns("genomic_region"), label = "Genomic region", value = NULL),
                   checkboxInput(ns("viewMode"), label = "Genomic View", value = FALSE),
                   checkboxInput(ns("useCustomRegions"), label = "Protein structures", value = FALSE),
                   checkboxInput(ns("other_tx"), label = "Full annotation", value = FALSE),
                   checkboxInput(ns("add_uorfs"), label = "uORF annotation", value = FALSE),
                   checkboxInput(ns("expression_plot"), label = "Add expression plot", value = FALSE),
                   checkboxInput(ns("summary_track"), label = "Summary top track", value = FALSE),
                   frame_type_select(ns, "summary_track_type", "Select summary display type"),
                   uiOutput(ns("clip")),
                   export_format_of_plot(ns)
          )
        ),
         actionButton(ns("go"), "Plot", icon = icon("rocket")), width=3
      )),
      mainPanel(
        jqui_resizable(plotlyOutput(outputId = ns("c"), height = "500px")) %>% shinycssloaders::withSpinner(color="#0dc5c1"),
        uiOutput(ns("variableUi"),
      ), plotlyOutput(outputId = ns("d")) %>% shinycssloaders::withSpinner(color="#0dc5c1"), width=9)
    )
  )
}

browser_server <- function(id, all_experiments, env, df, experiments,
                           tx, cds, libs, org, gene_name_list, rv,
                           browser_options) {
  moduleServer(
    id,
    function(input, output, session, all_exp = all_experiments) {
      study_and_gene_observers(input, output, session)
      output$clip <- renderUI({clipboard_url_button(input, session)})

      # Main plot controller, this code is only run if 'plot' is pressed
      mainPlotControls <- eventReactive(input$go,
        click_plot_browser_main_controller(input, tx, cds, libs, df),
        ignoreInit = check_plot_on_start(browser_options),
        ignoreNULL = FALSE)
      output$c <- renderPlotly(click_plot_browser(mainPlotControls, session)) %>%
        bindCache(mainPlotControls()$hash_browser) %>%
        bindEvent(mainPlotControls(), ignoreInit = FALSE, ignoreNULL = TRUE)

      # Protein display
      module_protein(input, output, gene_name_list, session)

      output$d <- renderPlotly({
        req(input$expression_plot == TRUE)
        click_plot_boxplot(mainPlotControls, session)}) %>%
          bindCache(mainPlotControls()$hash_expression)
      return(rv)
    }
  )
}


