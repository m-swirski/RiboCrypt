browser_allsamp_ui = function(id,  all_exp, browser_options,
                      metadata, gene_names_init, label = "Browser_allsamp") {
  ns <- NS(id)
  genomes <- unique(all_exp$organism)
  experiments <- all_exp$name
  normalizations <- normalizations("metabrowser")
  enrichment_test_types <- c(`Clusters (Order factor 1)` = "Clusters", `Ratio bins` = "Ratio bins", `Other gene tpm bins` = "Other gene tpm bins")
  columns_to_show <- c("BioProject", "CONDITION", "INHIBITOR",
                       "BATCH", "TIMEPOINT", "TISSUE", "CELL_LINE", "GENE", "FRACTION")
  metadata <- metadata[, colnames(metadata) %in% columns_to_show, with = FALSE]
  full_annotation <- as.logical(browser_options["full_annotation"])
  if (!is.null(browser_options["default_gene_meta"])) {
    all_isoforms <- subset(gene_names_init, label == browser_options["default_gene_meta"])
  } else all_isoforms <- NULL

  viewMode <- browser_options["default_view_mode"] == "genomic"
  introns_width <- as.numeric(browser_options["collapsed_introns_width"])
  panel_hidden_or_not_class <- ifelse(browser_options["hide_settings"] == "TRUE",
                                      "floating_settings_panel hidden",
                                      "floating_settings_panel")
  tabPanel(
    title = "MegaBrowser", icon = icon("chart-line"),
    # ---- HEAD with floating settings style ----
    browser_ui_settings_style(),
    # ---- Floating Settings Panel ----
    fluidRow(
      column(1, div(style = "position: relative;",
                    actionButton(ns("toggle_settings"), "", icon = icon("sliders-h"),
                                 style = "color: #fff; background-color: rgba(0,123,255,0.6); border-color: rgba(0,123,255,1); font-weight: bold;"),
                    # Floating settings panel overlays and drops down
                    div(id = ns("floating_settings"),
                        class = panel_hidden_or_not_class,  # Add a custom class
                        style = "position: absolute; top: 100%; left: 0; z-index: 10; background-color: white; padding: 10px; border: 1px solid #ddd; border-radius: 6px; width: max-content; min-width: 300px; box-shadow: 0px 4px 10px rgba(0,0,0,0.1);",
                        tabsetPanel(
                          tabPanel("MegaBrowser",
                                   fluidRow(
                                     organism_input_select(c("ALL", genomes), ns),
                                     experiment_input_select(experiments, ns, browser_options,
                                                             "default_experiment_meta"),
                                   ),
                                   fluidRow(column(6, metadata_input_select(ns, metadata, TRUE)),
                                            column(6, metadata_input_select(ns, metadata, FALSE,
                                                                            label = "Enrichment test on:",
                                                                            id = "enrichment_term",
                                                                            selected = enrichment_test_types[1],
                                                                            add = enrichment_test_types))),
                                   fluidRow(column(6, region_view_select(ns, "region_type", "Select region to view")),
                                            column(6, radioButtons(ns("plotType"), "Choose Plot Type:",
                                                choices = c("Plotly" = "plotly", "ComplexHeatmap" = "ggplot2"),
                                                selected = "ggplot2", inline = TRUE))),
                                   helper_button_redirect_call()
                          ),
                          tabPanel("Settings",
                                   fluidRow(
                                     numericInput(ns("extendLeaders"), "5' extension", 0),
                                     numericInput(ns("extendTrailers"), "3' extension", 0),
                                     numericInput(ns("collapsed_introns_width"), "Collapse Introns (nt flanks)",
                                                  introns_width)
                                   ),
                                   normalization_input_select(ns, choices = normalizations,
                                                              help_link = "mbrowser"),
                                   fluidRow(column(6, heatmap_color_select(ns)),
                                            column(6, sliderInput(ns("color_mult"), "Color scale zoom", min = 1, max = 10,
                                                                  value = 3))),
                                   fluidRow(column(6, gene_input_select(ns, label = "Sort by other gene", id = "other_gene")),
                                            column(6, textInput(ns("ratio_interval"), label = "Sort on interval/ratio : a:b;x:y", value = NULL))),
                                   checkboxInput(ns("display_annot"), label = "Display annotation", value = TRUE),
                                   checkboxInput(ns("summary_track"), label = "Summary top track", value = FALSE),
                                   checkboxInput(ns("frame"), label = "Split by frame", value = FALSE),
                                   frame_type_select(ns, "summary_track_type", "Select summary display type"),
                                   fluidRow(column(6, numericInput(ns("min_count"), "Minimum counts", min = 0, value = 100)),
                                            column(6, sliderInput(ns("kmer"), "K-mer length", min = 1, max = 20,
                                                                  value = as.numeric(browser_options["default_kmer"]))))
                          )
                        )
                    ))),
      column(2, gene_input_select(ns, FALSE, browser_options)),
      column(2, tx_input_select(ns, FALSE, all_isoforms, browser_options["default_isoform_meta"])),
      column(1, NULL, plot_button(ns("go"))),
      column(2, prettySwitch(ns("viewMode"), "Genomic View", value = viewMode,
                             status = "success", fill = TRUE, bigger = TRUE),
             prettySwitch(ns("other_tx"), "Full annotation", value = full_annotation,
                          status = "success", fill = TRUE, bigger = TRUE),
             prettySwitch(ns("collapsed_introns"), "Collapse introns", value = FALSE,
                          status = "success", fill = TRUE, bigger = TRUE)),
      column(2, motif_input_select(ns)),
      column(2, sliderInput(ns("clusters"), "K-means clusters", min = 1, max = 20,
                  value = 5))
    ),
    tags$hr(),
    # ---- Full Width Main Panel ----
    fluidRow(
      tabsetPanel(type = "tabs",
                  tabPanel("Heatmap", fluidRow(
                    jqui_resizable(
                      splitLayout(cellWidths = c("10%", "90%"),
                                  plotlyOutput(outputId = ns("d"), height = "618px", width = "130px"),
                                  uiOutput(outputId = ns("c")) %>%
                                    shinycssloaders::withSpinner(color="#0dc5c1"),
                                  width=9, cellArgs = list(style = "padding: 0px")))
                  ),
                  plotlyOutput(outputId = ns("e"))),
                  tabPanel("Statistics", DTOutput(outputId = ns("stats")) %>% shinycssloaders::withSpinner(color="#0dc5c1")),
                  tabPanel("Result table", DTOutput(outputId = ns("result_table")) %>% shinycssloaders::withSpinner(color="#0dc5c1"))
      )
    )
  )
}

browser_allsamp_server <- function(id, all_experiments, df, metadata,
                                   names_init, browser_options, exp_init,
                                   experiments = all_experiments$name) {
  moduleServer(
    id,
    function(input, output, session, all_exp = all_experiments) {
      ns <- NS(id)
      allsamples_observer_controller(input, output, session)
      plot_type <- "plotly"
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
                                                   isolate(input$color_mult),
                                                   isolate(input$plotType))) %>%
        bindEvent(table(), ignoreInit = FALSE, ignoreNULL = TRUE)

      output$myPlotlyPlot <- renderPlotly({
        req(input$plotType == "plotly")
        get_meta_browser_plot_full(table()$table,
                                   plot_object(), controller()$id,
                                   controller()$dff, controller()$summary_track,
                                   controller()$display_annot,
                                   region_type = controller()$region_type,
                                   isolate(input$plotType), controller()$tx_annotation,
                                   controller()$display_region,
                                   controller()$annotation,
                                   controller()$viewMode,
                                   controller()$collapsed_introns_width
                                   )}) %>%
        bindCache(controller()$table_plot) %>%
        bindEvent(plot_object(), ignoreInit = FALSE, ignoreNULL = TRUE)
      output$myGgplot <- renderPlot({
        req(input$plotType == "ggplot2")
        get_meta_browser_plot_full(table()$table,
                                   plot_object(), controller()$id,
                                   controller()$dff, controller()$summary_track,
                                   controller()$display_annot,
                                   region_type = controller()$region_type,
                                   isolate(input$plotType), controller()$tx_annotation,
                                   controller()$display_region,
                                   controller()$annotation,
                                   controller()$viewMode,
                                   controller()$collapsed_introns_width
                                   )}) %>%
        bindCache(controller()$table_plot) %>%
        bindEvent(plot_object(), ignoreInit = FALSE, ignoreNULL = TRUE)

      output$c <- renderUI(
        if (input$plotType == "plotly") {
          plotlyOutput(
            outputId = ns("myPlotlyPlot"), # Specific ID for Plotly output
            height = "700px",
            width = "100%"
          ) %>% shinycssloaders::withSpinner(color="#0dc5c1")
      } else if (input$plotType == "ggplot2") {
          plotOutput(
            outputId = ns("myGgplot"), # Specific ID for ggplot2 output
            height = "700px",
            width = "100%"
          ) %>% shinycssloaders::withSpinner(color="#0dc5c1")
      }) %>%
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
      output$e <- renderPlotly(allsamples_enrich_bar_plotly(meta_and_clusters()$enrich_dt)) %>%
        bindCache(controller()$table_hash, controller()$enrichment_term) %>%
        bindEvent(meta_and_clusters(),
                  ignoreInit = FALSE,
                  ignoreNULL = TRUE)
      output$stats <- renderDT(allsamples_meta_stats_shiny(meta_and_clusters()$enrich_dt)) %>%
        bindEvent(meta_and_clusters(),
                  ignoreInit = FALSE,
                  ignoreNULL = TRUE)
      output$result_table <- renderDT(cbind(attr(meta_and_clusters()$meta, "runIDs"),
                                            meta_and_clusters()$meta)[, c("index", "order") := NULL],
                                      extensions = 'Buttons',
                                      filter = "top",
                                      options = list(dom = 'Bfrtip',
                                                     buttons = NULL,
                                                     pageLength = 130)) %>%
        bindEvent(meta_and_clusters(),
                  ignoreInit = FALSE,
                  ignoreNULL = TRUE)
      observeEvent(input$toggle_settings, {
        # Toggle visibility by adding/removing 'hidden' class
        shinyjs::toggleClass(id = "floating_settings", class = "hidden")
      })
    }
  )
}
