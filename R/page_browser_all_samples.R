browser_allsamp_ui = function(id,  all_exp, browser_options,
                      metadata, label = "Browser_allsamp") {
  ns <- NS(id)
  genomes <- unique(all_exp$organism)
  experiments <- all_exp$name
  normalizations <- normalizations("metabrowser")
  tabPanel(
    title = "MetaBrowser", icon = icon("chart-line"),
    sidebarLayout(
      jqui_resizable(sidebarPanel(
        tabsetPanel(
          tabPanel("Browser",
                   organism_input_select(c("ALL", genomes), ns),
                   experiment_input_select(experiments, ns, browser_options,
                                           "default_experiment_meta"),
                   gene_input_select(ns, FALSE, browser_options),
                   metadata_input_select(ns, metadata = metadata),
                   sliderInput(ns("clusters"), "K-means clusters", min = 1, max = 5,
                               value = 1),
                   helper_button_redirect_call()
          ),
          tabPanel("Settings",
                   normalization_input_select(ns, choices = normalizations,
                                              help_link = "mbrowser"),
                   heatmap_color_select(ns),
                   sliderInput(ns("color_mult"), "Color scale zoom", min = 1, max = 10,
                               value = 3),
                   checkboxInput(ns("summary_track"), label = "Summary top track", value = FALSE),
                   frame_type_select(ns, "summary_track_type", "Select summary display type"),
                   tx_input_select(ns, FALSE),
                   sliderInput(ns("kmer"), "K-mer length", min = 1, max = 20,
                               value = as.numeric(browser_options["default_kmer"]))
          )
        ),
        actionButton(ns("go"), "Plot", icon = icon("rocket")), width=3
      )),
      mainPanel(
        fluidRow(
          splitLayout(cellWidths = c("10%", "90%"),
                      jqui_resizable(plotlyOutput(outputId = ns("d"), height = "700px", width = "130px")) %>%
                        shinycssloaders::withSpinner(color="#0dc5c1"),
                      jqui_resizable(plotOutput(outputId = ns("c"), height = "700px")) %>%
                        shinycssloaders::withSpinner(color="#0dc5c1"),
                      width=9)
        )
      )
    )
  )
}

browser_allsamp_server <- function(id, all_experiments, df, experiments,
                           metadata) {
  moduleServer(
    id,
    function(input, output, session, all_exp = all_experiments) {
      rv <- reactiveValues(lstval=isolate(df())@txdb,
                           curval=isolate(df())@txdb,
                           genome = "ALL",
                           exp = name(isolate(df())),
                           changed=FALSE)
      uses_libs <- FALSE
      org <- reactive("ALL")
      gene_name_list <- reactive(get_gene_name_categories_collection(df())) %>%
        bindCache(rv$curval) %>%
        bindEvent(rv$changed, ignoreNULL = TRUE)
      study_and_gene_observers(input, output, session)
      # Main plot controller, this code is only run if 'plot' is pressed
      mainPlotControls <- eventReactive(input$go,
                                        click_plot_browser_allsamp_controller(input, df),
                                        ignoreInit = TRUE,
                                        ignoreNULL = FALSE)
      # Main plot, this code is only run if 'plot' is pressed
      table <- reactive(compute_collection_table_shiny(mainPlotControls,
                                                   metadata = metadata)) %>%
        bindCache(mainPlotControls()$table_hash) %>%
        bindEvent(mainPlotControls()$table_hash,
                  ignoreInit = FALSE,
                  ignoreNULL = TRUE)

      plot_object <- reactive(get_meta_browser_plot(table(),
                                                   isolate(input$heatmap_color),
                                                   isolate(input$clusters),
                                                   isolate(input$color_mult))) %>%
        bindEvent(table(),
                  ignoreInit = FALSE,
                  ignoreNULL = TRUE)

      output$c <- renderPlot(plot_object()) %>%
        bindCache(mainPlotControls()$table_hash, input$heatmap_color,
                  isolate(input$color_mult)) %>%
        bindEvent(plot_object(),
                  ignoreInit = FALSE,
                  ignoreNULL = TRUE)
      output$d <- renderPlotly(allsamples_sidebar(mainPlotControls,
                                                  plot_object(),
                                                  metadata = metadata)) %>%
        bindCache(mainPlotControls()$table_hash) %>%
        bindEvent(plot_object(),
                  ignoreInit = FALSE,
                  ignoreNULL = TRUE)
    }
  )
}
