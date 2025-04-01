browser_ui = function(id,  all_exp, browser_options, gene_names_init,
                      libs, label = "Browser") {

  ns <- NS(id)
  genomes <- unique(all_exp$organism)
  experiments <- all_exp$name
  # TODO: move init_tx and init_libs to setup module
  all_isoforms <- subset(gene_names_init, label == browser_options["default_gene"])
  init_libs <- unlist(strsplit(browser_options["default_libs"], "\\|"))
  viewMode <- browser_options["default_view_mode"] == "genomic"
  copy_button_formatting <- tags$head(
    tags$style(HTML('#clip{background-color:orange}'))
  )

  tabPanel(
    tags$head(includeHTML(system.file("google_analytics_html",
                                      "google_analytics.html", package = "RiboCrypt"))),
    title = "browser", icon = icon("chart-line"),
    sidebarLayout(
      jqui_resizable(sidebarPanel(
        tabsetPanel(
          tabPanel("Browser",
                   fluidRow(column(6, organism_input_select(c("ALL", genomes), ns)),
                            column(6, experiment_input_select(experiments, ns, browser_options))),
                   fluidRow(column(6, gene_input_select(ns, FALSE, browser_options)),
                           column(6, tx_input_select(ns, FALSE, all_isoforms, browser_options["default_isoform"]))),
                   library_input_select(ns, TRUE, libs, init_libs),
                   fluidRow(column(6, frame_type_select(ns, selected =
                                       browser_options["default_frame_type"])),
                            column(6, sliderInput(ns("kmer"), "K-mer length", min = 1, max = 20,
                          value = as.numeric(browser_options["default_kmer"])))),
                   shinyjs::useShinyjs(),
                   rclipboardSetup(),
                   copy_button_formatting
          ),
          tabPanel("Settings",
                   fluidRow(column(6, numericInput(ns("extendLeaders"), "5' extension", 0)),
                            column(6, numericInput(ns("extendTrailers"), "3' extension", 0))),
                   textInput(ns("customSequence"), label = "Custom sequences highlight", value = NULL),
                   textInput(ns("genomic_region"), label = "Genomic region", value = NULL),
                   textInput(ns("zoom_range"), label = "Zoom interval", value = NULL),
                   fluidRow(column(4, checkboxInput(ns("other_tx"), label = "Full annotation", value = FALSE)),
                            column(4, checkboxInput(ns("add_uorfs"), label = "uORF annotation", value = FALSE)),
                            column(4, checkboxInput(ns("add_translon"), label = "Predicted translons", value = FALSE))),
                   checkboxInput(ns("log_scale"), label = "Log Scale", value = FALSE),
                   fluidRow(column(4, checkboxInput(ns("expression_plot"), label = "Gene expression plot", value = FALSE)),
                            column(4, checkboxInput(ns("useCustomRegions"), label = "Protein structures", value = TRUE)),
                            column(4, checkboxInput(ns("phyloP"), label = "Conservation (phyloP)", value = FALSE))),
                   fluidRow(column(6, checkboxInput(ns("withFrames"), label = "Split color Frames", value = TRUE)),
                            column(6, frame_subsetter_select(ns))),
                   fluidRow(column(6, checkboxInput(ns("summary_track"), label = "Summary top track", value = FALSE)),
                            column(6, frame_type_select(ns, "summary_track_type", "Select summary display type"))),
                   uiOutput(ns("clip")),
                   export_format_of_plot(ns)
          )
        ),
        tags$hr(style = "border: 1px solid black; margin-top: 0px; margin-bottom: 10px;"),
        fluidRow(column(8, actionButton(ns("go"), "Generate Plot", icon = icon("rocket"),
                                        class = "btn-primary", style = "width: 100%; font-size: 16px;")),
                 column(4, prettySwitch(inputId = ns("viewMode"), label = "Genomic View", value = viewMode,
                   status = "success", fill = TRUE, bigger = TRUE))
        ), width=3
      )),
      mainPanel(
        jqui_resizable(plotlyOutput(outputId = ns("c"), height = "500px")) %>% shinycssloaders::withSpinner(color="#0dc5c1"),
        plotlyOutput(outputId = ns("e"), height = "50px"),
        uiOutput(ns("variableUi"),),
        plotlyOutput(outputId = ns("d")) %>% shinycssloaders::withSpinner(color="#0dc5c1"), width=9)
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

      bottom_panel <- reactive(bottom_panel_shiny(mainPlotControls))  %>%
        bindCache(mainPlotControls()$hash_bottom) %>%
        bindEvent(mainPlotControls(), ignoreInit = FALSE, ignoreNULL = TRUE)

      output$c <- renderPlotly(browser_track_panel_shiny(mainPlotControls, bottom_panel(), session)) %>%
        bindCache(mainPlotControls()$hash_browser) %>%
        bindEvent(bottom_panel(), ignoreInit = FALSE, ignoreNULL = TRUE)


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


