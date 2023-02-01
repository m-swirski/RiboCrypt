browser_ui = function(id, label = "Browser", all_exp) {
  ns <- NS(id)
  genomes <- unique(all_exp$organism)
  experiments <- all_exp$name
  tabPanel(
    title = "browser", icon = icon("chart-line"),
    sidebarLayout(
      jqui_resizable(jqui_draggable(sidebarPanel(
        tabsetPanel(
          tabPanel("Browser",
                   organism_input_select(c("ALL", genomes), ns),
                   experiment_input_select(experiments, ns),
                   gene_input_select(ns),
                   tx_input_select(ns),
                   library_input_select(ns),
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
                   checkboxInput(ns("summary_track"), label = "Summary top track", value = FALSE),
                   frame_type_select(ns, "summary_track_type", "Select summary display type"),
                   export_format_of_plot(ns)
          ),
        ),
         actionButton(ns("go"), "Plot", icon = icon("rocket")),
      ))),
      mainPanel(
        jqui_resizable(plotlyOutput(outputId = ns("c"), height = "500px")) %>% shinycssloaders::withSpinner(color="#0dc5c1"),
        uiOutput(ns("variableUi"))
      )
    )
  )
}

browser_server <- function(id, all_experiments, env) {
  moduleServer(
    id,
    function(input, output, session, all_exp = all_experiments) {
      # Organism / study objects
      org_and_study_changed_checker(input, output, session)
      # Gene objects
      tx <- reactive({loadRegion(df())}) %>%
        bindEvent(rv_changed(), ignoreNULL = TRUE)
      cds <- reactive(loadRegion(df(), part = "cds")) %>%
        bindEvent(rv_changed(), ignoreNULL = TRUE)
      gene_name_list <- reactive(get_gene_name_categories(df())) %>%
        bindCache(rv$curval) %>%
        bindEvent(rv_changed(), ignoreNULL = TRUE)
      # Update main side panels
      study_and_gene_observers(input, output, session)

      # Main plot controller, this code is only run if 'plot' is pressed
      mainPlotControls <- eventReactive(input$go,
        click_plot_browser_main_controller(input, tx, cds, libs, df))
      # Main plot, this code is only run if 'plot' is pressed
      output$c <- renderPlotly(click_plot_browser(mainPlotControls, session))
      # Protein display
      module_protein(input, output, session)
    }
  )
}
