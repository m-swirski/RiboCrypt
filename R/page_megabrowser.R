browser_allsamp_ui = function(id,  all_exp, browser_options,
                      metadata, label = "Browser_allsamp") {
  ns <- NS(id)
  genomes <- unique(all_exp$organism)
  experiments <- all_exp$name
  normalizations <- normalizations("metabrowser")
  enrichment_test_types <- c("Clusters", "Ratio bins", "Other gene tpm bins")
  columns_to_show <- c("BioProject", "CONDITION", "INHIBITOR",
                       "BATCH", "TIMEPOINT", "TISSUE", "CELL_LINE", "GENE", "FRACTION")
  metadata <- metadata[, colnames(metadata) %in% columns_to_show, with = FALSE]
  tabPanel(
    title = "MegaBrowser", icon = icon("chart-line"),
    sidebarLayout(
      jqui_resizable(sidebarPanel(
        tabsetPanel(
          tabPanel("Browser",
                   organism_input_select(c("ALL", genomes), ns),
                   experiment_input_select(experiments, ns, browser_options,
                                           "default_experiment_meta"),
                   gene_input_select(ns, FALSE, browser_options),
                   metadata_input_select(ns, metadata, TRUE),
                   metadata_input_select(ns, metadata, FALSE,
                                         label = "Enrichment test on:",
                                         id = "enrichment_term",
                                         add = enrichment_test_types),
                   sliderInput(ns("clusters"), "K-means clusters", min = 1, max = 20,
                               value = 5),
                   helper_button_redirect_call()
          ),
          tabPanel("Settings",
                   normalization_input_select(ns, choices = normalizations,
                                              help_link = "mbrowser"),
                   heatmap_color_select(ns),
                   sliderInput(ns("color_mult"), "Color scale zoom", min = 1, max = 10,
                               value = 3),
                   gene_input_select(ns, label = "Sort by other gene", id = "other_gene"),
                   motif_input_select(ns),
                   textInput(ns("ratio_interval"), label = "Sort on interval/ratio : a:b;x:y", value = NULL),
                   checkboxInput(ns("display_annot"), label = "Display annotation", value = TRUE),
                   checkboxInput(ns("summary_track"), label = "Summary top track", value = FALSE),
                   checkboxInput(ns("frame"), label = "Split by frame", value = FALSE),
                   frame_type_select(ns, "summary_track_type", "Select summary display type"),
                   region_view_select(ns, "region_type", "Select region to view"),
                   tx_input_select(ns, FALSE),
                   numericInput(ns("min_count"), "Minimum counts", min = 0, value = 100),
                   sliderInput(ns("kmer"), "K-mer length", min = 1, max = 20,
                               value = as.numeric(browser_options["default_kmer"]))
          )
        ),
        actionButton(ns("go"), "Plot", icon = icon("rocket")), width=3
      )),
      mainPanel(
        tabsetPanel(type = "tabs",
          tabPanel("Heatmap", fluidRow(
            jqui_resizable(
              splitLayout(cellWidths = c("10%", "90%"),
                        plotlyOutput(outputId = ns("d"), height = "618px", width = "130px"),
                        plotOutput(outputId = ns("c"), height = "700px", width = "105%") %>%
                          shinycssloaders::withSpinner(color="#0dc5c1"),
                        width=9, cellArgs = list(style = "padding: 0px")))
        ),
        plotlyOutput(outputId = ns("e"))),
        tabPanel("Statistics", DTOutput(outputId = ns("stats")) %>% shinycssloaders::withSpinner(color="#0dc5c1"))
      ))
    )
  )
}

browser_allsamp_server <- function(id, all_experiments, df, metadata,
                                   names_init, browser_options,
                                   experiments = all_experiments$name) {
  moduleServer(
    id,
    function(input, output, session, all_exp = all_experiments) {
      allsamples_observer_controller(input, output, session)

      # Main plot controller, this code is only run if 'plot' is pressed
      controller <- eventReactive(input$go,
                                  click_plot_browser_allsamp_controller(input, df, gene_name_list),
                                  ignoreInit = TRUE,
                                  ignoreNULL = FALSE)
      # Main plot, this code is only run if 'plot' is pressed
      table <- reactive(compute_collection_table_shiny(controller,
                                                   metadata = metadata)) %>%
        bindCache(controller()$table_hash) %>%
        bindEvent(controller()$table_hash, ignoreInit = FALSE, ignoreNULL = TRUE)

      plot_object <- reactive(get_meta_browser_plot(table()$table,
                                                   isolate(input$heatmap_color),
                                                   isolate(input$clusters),
                                                   isolate(input$color_mult))) %>%
        bindEvent(table(), ignoreInit = FALSE, ignoreNULL = TRUE)

      output$c <- renderPlot(get_meta_browser_plot_full(table()$table,
                        plot_object(), controller()$id,
                        controller()$dff, controller()$summary_track,
                        controller()$display_annot,
                        region_type = controller()$region_type
                             )) %>%
        bindCache(controller()$table_hash, input$heatmap_color,
                  isolate(input$color_mult)) %>%
        bindEvent(plot_object(), ignoreInit = FALSE, ignoreNULL = TRUE)

      # Addition plots
      meta_and_clusters <- reactive(
          allsamples_metadata_clustering(table()$metadata_field, plot_object(),
                                         controller()$enrichment_term)) %>%
        bindCache(controller()$table_hash, controller()$enrichment_term) %>%
        bindEvent(plot_object(), ignoreInit = FALSE, ignoreNULL = TRUE)

      output$d <- renderPlotly(allsamples_sidebar(meta_and_clusters()$meta)) %>%
        bindCache(controller()$table_hash) %>%
        bindEvent(meta_and_clusters(),
                  ignoreInit = FALSE,
                  ignoreNULL = TRUE)
      output$e <- renderPlotly(allsamples_enrich_bar_plot(meta_and_clusters()$enrich_dt)) %>%
        bindCache(controller()$table_hash, controller()$enrichment_term) %>%
        bindEvent(meta_and_clusters(),
                  ignoreInit = FALSE,
                  ignoreNULL = TRUE)
      output$stats <- renderDT(allsamples_meta_stats_shiny(meta_and_clusters()$enrich_dt)) %>%
        bindEvent(meta_and_clusters(),
                  ignoreInit = FALSE,
                  ignoreNULL = TRUE)
    }
  )
}
