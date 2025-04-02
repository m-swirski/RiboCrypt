heatmap_ui <- function(id, all_exp, browser_options, libs, label = "Heatmap") {
  ns <- NS(id)
  genomes <- unique(all_exp$organism)
  experiments <- all_exp$name
  tabPanel(
    title = "Motif metaplot", icon = icon("layer-group"),
    sidebarLayout(
      sidebarPanel(
        tabsetPanel(
          tabPanel("Motif metaplot",
                   fluidRow(column(6, organism_input_select(c("ALL", genomes), ns)),
                            column(6, experiment_input_select(experiments, ns, browser_options))),
                   fluidRow(column(6, gene_input_select(ns)),
                            column(6, tx_input_select(ns))),
                   library_input_select(ns, FALSE, libs),
                   fluidRow(column(6, heatmap_region_select(ns, "User defined Motif")),
                            column(3, numericInput(ns("extendLeaders"), "5' extension", 30)),
                            column(3, numericInput(ns("extendTrailers"), "3' extension", 30))),
                   textInput(ns("customSequence"), label = "Anchor on motif", value = NULL),
                   normalization_input_select(ns),
                   fluidRow(column(6, numericInput(ns("readlength_min"), "Min Readlength", 26)),
                            column(6, numericInput(ns("readlength_max"), "Max Readlength", 34)))),
          tabPanel("Settings",
                   checkboxInput(ns("summary_track"), label = "Summary top track", value = FALSE),
                   checkboxInput(ns("p_shifted"), label = "p_shifted", value = TRUE))),
        plot_button(ns("go"))
      ),
      mainPanel(
        tabsetPanel(type = "tabs",
                    tabPanel("Heatmap", plotlyOutput(outputId = ns("c"), height = "500px") %>% shinycssloaders::withSpinner(color="#0dc5c1"),
                             uiOutput(ns("variableUi"))),
                    tabPanel("Shift table", tableOutput(outputId = ns("shift_table")) %>% shinycssloaders::withSpinner(color="#0dc5c1"))))
  ))
}

heatmap_server <- function(id, all_experiments, env, df, experiments, tx, cds,
                           libs, org, gene_name_list, rv) {
  moduleServer(
    id,
    function(input, output, session, all_exp = all_experiments) {
      # Update main side panels
      all_is_gene <- TRUE
      study_and_gene_observers(input, output, session)

      # Main plot, this code is only run if 'plot' is pressed
      mainPlotControls <- eventReactive(input$go,
        click_plot_heatmap_main_controller(input, tx, cds, libs, df))

      anchor_points <- reactive(anchor_points_shiny(mainPlotControls)) %>%
        bindCache(mainPlotControls()$hash_string_anchor) %>%
        bindEvent(mainPlotControls(), ignoreNULL = TRUE)

      coverage <- reactive(heatmap_data(mainPlotControls, tx, anchor_points())) %>%
        bindCache(mainPlotControls()$hash_string) %>%
        bindEvent(anchor_points(), ignoreNULL = TRUE)

      output$c <- renderPlotly(reactive_heatmap_plot(mainPlotControls, coverage())) %>%
        bindEvent(coverage(), ignoreNULL = TRUE)
      output$shift_table <- renderTable(mainPlotControls()$shift_table)
      return(rv)
  }
  )
}
