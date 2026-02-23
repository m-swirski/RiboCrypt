browser_allsamp_ui = function(id,  all_exp, browser_options,
                      metadata, gene_names_init, label = "Browser_allsamp") {
  ns <- NS(id)
  genomes <- unique(all_exp$organism)
  experiments <- all_exp$name
  normalizations <- normalizations("metabrowser")
  enrichment_test_types <- c(`Clusters (Order factor 1)` = "Clusters", `Ratio bins` = "Ratio bins", `Other gene tpm bins` = "Other gene tpm bins")
  columns_to_show <- c("BioProject", "CONDITION", "INHIBITOR",
                       "BATCH", "TIMEPOINT", "TISSUE", "CELL_LINE", "GENE", "FRACTION")
  metadata <- metadata[, ..columns_to_show]
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
                                     column(6, organism_input_select(c("ALL", genomes), ns)),
                                     column(6, experiment_input_select(experiments, ns, browser_options,
                                                             "default_experiment_meta")),
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
                                                selected = "plotly", inline = TRUE))),
                                   tags$hr(style = "padding-top: 50px; padding-bottom: 50px;"),
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
      tabsetPanel(id = ns("mb_tabs"), type = "tabs",
                  tabPanel("Heatmap", fluidRow(
                    jqui_resizable(
                      div(
                        style = paste(
                          "display: grid;",
                          "grid-template-columns: 8% 92%;",
                          "grid-template-rows: 15% 75% 10%;",
                          "height: 700px;",
                          "width: 100%;",
                          "gap: 0;"
                        ),
                        div(
                          style = "grid-column: 1; grid-row: 2; overflow: visible;",
                          plotly::plotlyOutput(outputId = ns("d"),
                                               height = "100%", width = "100%")
                        ),
                        div(
                          style = "grid-column: 2; grid-row: 1 / span 3;",
                          uiOutput(outputId = ns("c")) %>%
                            shinycssloaders::withSpinner(color = "#0dc5c1")
                        )
                      )
                    )
                  ),
                  plotlyOutput(outputId = ns("e"))),
                  tabPanel("Statistics", DTOutput(outputId = ns("stats")) %>% shinycssloaders::withSpinner(color="#0dc5c1")),
                  tabPanel("Result table",
                           uiOutput(outputId = ns("result_table_controls")),
                           DTOutput(outputId = ns("result_table")) %>%
                             shinycssloaders::withSpinner(color="#0dc5c1"))
      )
    )
  )
}

browser_allsamp_server <- function(id, all_experiments, df, metadata,
                                   names_init, browser_options, exp_init,
                                   exps_dir, gg_theme, experiments = all_experiments$name) {
  moduleServer(
    id,
    function(input, output, session, all_exp = all_experiments) {
      ns <- NS(id)
      allsamples_observer_controller(input, output, session)
      plot_type <- "plotly"
      # Main plot controller, this code is only run if 'plot' is pressed
      controller <- reactive(mb_controller_shiny(input, df, gene_name_list)) %>%
        bindCache(input_to_list(input)) %>%
        bindEvent(input$go, ignoreInit = TRUE, ignoreNULL = FALSE)
      # Table
      table <- reactive(compute_collection_table_shiny(controller, metadata = metadata)) %>%
        bindCache(controller()$table_hash) %>%
        bindEvent(controller()$table_hash, ignoreInit = FALSE, ignoreNULL = TRUE)

      # Heatmap (middle right)
      plot_object <- reactive(mb_plot_object_shiny(table()$table, input)) %>%
        bindCache(controller()$table_plot_hash) %>%
        bindEvent(table(), ignoreInit = FALSE, ignoreNULL = TRUE)

      mb_top_plot <- reactive(mb_top_plot_shiny(table()$table)) %>%
        bindCache(controller()$table_hash) %>%
        bindEvent(plot_object(), ignoreInit = FALSE, ignoreNULL = TRUE)

      mb_mid_plot <- reactive(mb_mid_plot_shiny(plot_object(), input$plotType)) %>%
        bindCache(controller()$table_plot_hash) %>%
        bindEvent(plot_object(), ignoreInit = FALSE, ignoreNULL = TRUE)

      mb_bottom_plot <- reactive(mb_bottom_plot_shiny(controller, gg_theme)) %>%
        bindCache(controller()$table_hash) %>%
        bindEvent(plot_object(), ignoreInit = FALSE, ignoreNULL = TRUE)

      output$myPlotlyPlot <- renderPlotly({
        req(input$plotType == "plotly")
        mb_mid_plot()
      }) %>%
        bindCache(controller()$table_plot_hash) %>%
        bindEvent(plot_object(), ignoreInit = FALSE, ignoreNULL = TRUE)

      mb_mid_image <- reactive(mb_mid_image_shiny(plot_object(), session, ns)) %>%
        bindCache(controller()$table_plot_hash) %>%
        bindEvent(plot_object(), ignoreInit = FALSE, ignoreNULL = TRUE)

      output$mb_subplot_static <- renderPlotly({
        req(input$plotType == "ggplot2")
        mb_static_subplot_shiny(mb_top_plot(), mb_mid_image(), mb_bottom_plot())
      }) %>%
        bindCache(controller()$table_plot_hash) %>%
        bindEvent(plot_object(), ignoreInit = FALSE, ignoreNULL = TRUE)

      output$mb_top_summary <- renderPlotly({
        mb_top_plot()
      }) %>%
        bindCache(controller()$table_hash) %>%
        bindEvent(plot_object(), ignoreInit = FALSE, ignoreNULL = TRUE)

      output$mb_bottom_gene <- renderPlotly({
        mb_bottom_plot()
      }) %>%
        bindCache(controller()$table_hash) %>%
        bindEvent(plot_object(), ignoreInit = FALSE, ignoreNULL = TRUE)

      observe({
        req(input$plotType == "plotly")
        ed <- suppressWarnings(plotly::event_data("plotly_relayout", source = "mb_mid"))
        req(!is.null(ed))
        y_max <- ncol(table()$table)
        sync_megabrowser_x_shiny(ed, session, y_max = y_max, y_reversed = TRUE)
      })

      output$c <- renderUI(renderMegabrowser(input$plotType, ns)) %>%
        bindCache(controller()$table_plot_hash) %>%
        bindEvent(plot_object(), ignoreInit = FALSE, ignoreNULL = TRUE)

      # Additional plots and tables
      meta_and_clusters <- reactive(
        allsamples_metadata_clustering(table(),
                                       controller()$enrichment_term)) %>%
        bindCache(controller()$table_hash, controller()$enrichment_term) %>%
        bindEvent(plot_object(), ignoreInit = FALSE, ignoreNULL = TRUE)

      selected_enrich_filters <- reactiveVal(NULL)
      observeEvent(meta_and_clusters(), {
        selected_enrich_filters(NULL)
      }, ignoreInit = TRUE)

      output$d <- renderPlotly({
        allsamples_sidebar_plotly(meta_and_clusters()$meta)
      }) %>%
        bindCache(controller()$table_hash) %>%
        bindEvent(meta_and_clusters(),
                  ignoreInit = FALSE,
                  ignoreNULL = TRUE)
      output$e <- renderPlotly(allsamples_enrich_bar_plotly(meta_and_clusters()$enrich_dt)) %>%
        bindCache(controller()$table_hash, controller()$enrichment_term) %>%
        bindEvent(meta_and_clusters(),
                  ignoreInit = FALSE,
                  ignoreNULL = TRUE)

      observeEvent(plotly::event_data("plotly_click", source = "mb_enrich"), {
        ed <- plotly::event_data("plotly_click", source = "mb_enrich")
        req(!is.null(ed))

        category <- as.character(ed$x)
        cluster <- NULL
        if (!is.null(ed$customdata) && !identical(ed$customdata, "counts")) {
          cluster <- as.character(ed$customdata)
        }

        selected_enrich_filters(list(category = category, cluster = cluster))
        updateTabsetPanel(session, "mb_tabs", selected = "Result table")
        shinyjs::runjs(
          sprintf(
            "document.getElementById('%s').scrollIntoView({behavior:'smooth'});",
            ns("result_table")
          )
        )
      })

      selected_cluster_from_filtered <- reactive({
        filt <- selected_enrich_filters()
        if (is.null(filt)) return(NULL)
        if (!is.null(filt$cluster)) return(as.character(filt$cluster))

        tbl <- filtered_meta_table()
        if ("cluster" %in% names(tbl)) {
          vals <- unique(tbl$cluster)
        } else if ("Cluster" %in% names(tbl)) {
          vals <- unique(tbl$Cluster)
        } else {
          vals <- character(0)
        }
        if (length(vals) == 1) return(as.character(vals))
        NULL
      })

      output$result_table_controls <- renderUI({
        filt <- selected_enrich_filters()
        if (is.null(filt)) return(NULL)

        show_cluster_btn <- FALSE
        cluster_val <- selected_cluster_from_filtered()
        if (!is.null(cluster_val)) {
          already_full_cluster <- is.null(filt$category) && !is.null(filt$cluster)
          show_cluster_btn <- !already_full_cluster
        }

        tagList(
          actionButton(ns("reset_result_table"), "Show full table"),
          if (show_cluster_btn) actionButton(ns("expand_cluster"), "Show full cluster") else NULL
        )
      })

      observeEvent(input$reset_result_table, {
        selected_enrich_filters(NULL)
      }, ignoreInit = TRUE)

      observeEvent(input$expand_cluster, {
        cluster_val <- selected_cluster_from_filtered()
        req(!is.null(cluster_val))
        selected_enrich_filters(list(category = NULL, cluster = cluster_val))
      }, ignoreInit = TRUE)

      filtered_meta_table <- reactive({
        tbl <- allsamples_meta_table(meta_and_clusters())
        filt <- selected_enrich_filters()
        if (is.null(filt)) return(tbl)
        if (!is.null(filt$category) && "grouping" %in% names(tbl)) {
          tbl <- tbl[as.character(grouping) %in% filt$category]
        }
        if (!is.null(filt$cluster)) {
          if ("cluster" %in% names(tbl)) {
            tbl <- tbl[as.character(cluster) %in% filt$cluster]
          } else if ("Cluster" %in% names(tbl)) {
            tbl <- tbl[as.character(Cluster) %in% filt$cluster]
          }
        }
        tbl
      })

      output$stats <- renderDT(allsamples_meta_stats_shiny(meta_and_clusters()$enrich_dt)) %>%
        bindEvent(meta_and_clusters(),
                  ignoreInit = FALSE,
                  ignoreNULL = TRUE)

      output$result_table <- renderDT(filtered_meta_table(),
                                      extensions = 'Buttons', filter = "top",
                                      options = list(dom = 'Bfrtip',
                                                     buttons = NULL,
                                                     pageLength = 130)) %>%
        bindEvent(meta_and_clusters(), selected_enrich_filters(),
                  ignoreInit = FALSE,
                  ignoreNULL = TRUE)
      observeEvent(input$toggle_settings, {
        # Toggle visibility by adding/removing 'hidden' class
        shinyjs::toggleClass(id = "floating_settings", class = "hidden")
      })
    }
  )
}
