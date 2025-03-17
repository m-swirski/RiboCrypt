quality_ui <- function(id, all_exp, browser_options, libs, label = "Periodicity") {
  ns <- NS(id)
  genomes <- unique(all_exp$organism)
  experiments <- all_exp$name
  tabPanel(
    title = "Frame bias", icon = icon("layer-group"),
    sidebarLayout(
      sidebarPanel(
        tabsetPanel(
          tabPanel("Periodicity",
                   organism_input_select(c("ALL", genomes), ns),
                   experiment_input_select(experiments, ns, browser_options),
                   gene_input_select(ns),
                   tx_input_select(ns),
                   library_input_select(ns, FALSE, libs),
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

quality_server <- function(id, all_experiments, env, df, experiments,
                           tx, cds, libs, org, gene_name_list, rv) {
  moduleServer(
    id,
    function(input, output, session, all_exp = all_experiments) {
      # Update main side panels
      length_table <- reactive(optimizedTranscriptLengths(df(), TRUE, TRUE)) %>%
        bindCache(rv$curval) %>%
        bindEvent(rv$changed, ignoreNULL = TRUE)

      all_is_gene <- TRUE
      study_and_gene_observers(input, output, session)

      # Main plot, this code is only run if 'plot' is pressed
      mainPlotControls <- eventReactive(input$go,
        click_plot_heatmap_main_controller(input, tx, cds, libs, df,
                                           length_table))

      coverage <- reactive(heatmap_data(mainPlotControls, tx, length_table)) %>%
        bindCache(mainPlotControls()$extendLeaders, mainPlotControls()$extendTrailers,
                  mainPlotControls()$normalization, mainPlotControls()$region,
                  ORFik:::name_decider(mainPlotControls()$dff, naming = "full"),
                  mainPlotControls()$readlength_min,
                  mainPlotControls()$readlength_max)



      output$c <- renderPlotly({
        coverage()[, frame := as.factor(position %% 3), by = fraction]
        return(periodicity_plot(coverage(), fft = TRUE))
      })
      return(rv)
    }
  )
}
