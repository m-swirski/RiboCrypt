browser_allsamp_ui = function(id,  all_exp, browser_options,
                      metadata, gene_names_init, label = "Browser_allsamp") {
  ns <- NS(id)
  genomes <- unique(all_exp$organism)
  experiments <- all_exp$name
  normalizations <- normalizations("metabrowser")
  enrichment_test_types <- c(`Clusters (Order factor 1)` = "Clusters", `Ratio bins` = "Ratio bins", `Other gene tpm bins` = "Other gene tpm bins")
  if (!is.null(metadata)) {
    columns_to_show <- c("BioProject", "YEAR","CONDITION", "INHIBITOR",
                         "BATCH", "TIMEPOINT", "TISSUE", "CELL_LINE", "GENE", "FRACTION",
                         "Cancer_type", "Cell_model", "Cell_type")
    columns_to_show <- columns_to_show[columns_to_show %in% colnames(metadata)]
    metadata <- metadata[, ..columns_to_show]
  }

  full_annotation <- as.logical(browser_options["full_annotation"])
  translons <- isTRUE(as.logical(browser_options["translons"]))
  translons_transcode <- isTRUE(as.logical(browser_options["translons_transcode"]))
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
    megabrowser_layout_assets(),
    div(id = ns("workspace"), class = "mega-workspace",
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
                                   fluidRow(
                                     column(4, checkboxInput(ns("add_translon"), "Predicted translons (Our all-merged: T)", translons)),
                                     column(4, checkboxInput(ns("add_translons_transcode"), "Predicted translons (TransCode: TC)",
                                                             translons_transcode))
                                   ),
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
    megabrowser_view_controls(ns),
    # ---- Full Width Main Panel ----
    fluidRow(
      tabsetPanel(id = ns("mb_tabs"), type = "tabs",
                  tabPanel("Heatmap", fluidRow(
                    jqui_resizable(
                      div(
                        class = "mega-plot-grid",
                        div(
                          class = "mega-plot-sidebar",
                          plotly::plotlyOutput(outputId = ns("d"),
                                               height = "100%", width = "100%")
                        ),
                        div(
                          class = "mega-plot-main",
                          uiOutput(outputId = ns("c")) %>%
                            shinycssloaders::withSpinner(color = "#0dc5c1")
                        )
                      )
                    )
                  )),
                  tabPanel("Factor enrichment", tabsetPanel(id = ns("factor_tabs"),
                  tabPanel("Enrichment", plotlyOutput(outputId = ns("e"))),
                  tabPanel("Group summary",
                           div(class = "mega-layout-controls",
                               downloadButton(ns("download_groups"), "Group summary CSV"),
                               downloadButton(ns("download_membership"), "Library membership CSV")),
                           DTOutput(ns("group_summary"))),
                  tabPanel("Statistics",
                           DTOutput(outputId = ns("stats")) %>% shinycssloaders::withSpinner(color="#0dc5c1")),
                  tabPanel("Result table",
                           uiOutput(outputId = ns("result_table_controls")),
                           DTOutput(outputId = ns("result_table")) %>%
                             shinycssloaders::withSpinner(color="#0dc5c1")))),
                  megabrowser_translon_ui(ns)
      )
    ))
  )
}

#' Full x reset range for the active megabrowser state.
#' @noRd
megabrowser_reset_range_shiny <- function(controller, table) {
  megabrowser_full_x_range(controller()$display_region, table()$table)
}

#' Reset layout for the central megabrowser heatmap.
#' @noRd
megabrowser_mid_reset_layout <- function(full_range, table_obj) {
  list(
    "xaxis.range" = full_range,
    "xaxis.autorange" = FALSE,
    "yaxis.range" = c(0.5, ncol(table_obj$table) + 0.5),
    "yaxis.autorange" = FALSE
  )
}

#' Reset layout for top and bottom megabrowser peer tracks.
#' @noRd
megabrowser_peer_reset_layout <- function(full_range) {
  list("xaxis.range" = full_range, "xaxis.autorange" = FALSE)
}

#' Add double-click reset behavior to the central megabrowser heatmap.
#' @noRd
megabrowser_mid_reset_plot <- function(plot, controller, table, ns) {
  full_range <- megabrowser_reset_range_shiny(controller, table)
  addMegabrowserDoubleClickReset(
    plot,
    reset_range = full_range,
    peer_ids = c(ns("mb_top_summary"), ns("mb_bottom_gene")),
    reset_layout = megabrowser_mid_reset_layout(full_range, table()),
    peer_reset_layout = megabrowser_peer_reset_layout(full_range)
  )
}

#' Add double-click reset behavior to a megabrowser peer track.
#' @noRd
megabrowser_peer_reset_plot <- function(plot, controller, table, peer_ids) {
  full_range <- megabrowser_reset_range_shiny(controller, table)
  addMegabrowserDoubleClickReset(plot, reset_range = full_range, peer_ids = peer_ids)
}

browser_allsamp_server <- function(id, all_exp, df, experiments,
                                   gene_name_list, tx, cds, org, motif_name_list,
                                   metadata, browser_options, rv,
                                   templates = NULL
) {
  moduleServer(
    id,
    function(input, output, session) {
      ns <- NS(id)
      allsamples_observer_controller(input, output, session)
      plot_type <- "plotly"
      # Main plot controller, this code is only run if 'plot' is pressed
      controller <- reactive(mb_controller_shiny(input, df, gene_name_list, cds, tx)) %>%
        bindCache(browser_controller_cache_version(), input_to_list(input)) %>%
        bindEvent(input$go, ignoreInit = TRUE, ignoreNULL = FALSE)
      # Table
      table <- reactive(compute_collection_table_shiny(controller, metadata = metadata)) %>%
        bindCache(controller()$table_hash) %>%
        bindEvent(controller()$table_hash, ignoreInit = FALSE, ignoreNULL = TRUE)
      grouped_metadata <- reactive(allsamples_metadata_clustering(table(), controller()$enrichment_term, compute_stats = FALSE)) %>%
        bindCache(controller()$table_hash, controller()$enrichment_term)
      selected_groups <- reactive({
        available <- names(megabrowser_display_groups(grouped_metadata()$meta, colnames(table()$table)))
        selected <- intersect(available, input$visible_groups)
        if (!length(selected) || length(selected) == length(available)) character() else selected
      })
      observeEvent(table(), {
        groups <- megabrowser_display_groups(grouped_metadata()$meta, colnames(table()$table))
        labels <- paste0(names(groups), " (", format(lengths(groups), big.mark = ","), " libraries)")
        updateSelectizeInput(session, "visible_groups", choices = stats::setNames(names(groups), labels),
                             selected = character(), server = TRUE)
      })
      display_table <- reactive({
        focused <- megabrowser_focused_display(table()$table, grouped_metadata()$meta, selected_groups())
        if (isTRUE(input$collapsed_clusters)) megabrowser_collapsed_display(focused$table, focused$meta) else focused
      }) %>% bindCache(controller()$table_hash, controller()$enrichment_term, isTRUE(input$collapsed_clusters), selected_groups())
      output$display_counts <- renderText(megabrowser_display_counts(table()$table, display_table()))
      megabrowser_group_outputs(output, table, grouped_metadata, selected_groups)
      megabrowser_translon_outputs(input, output, session, controller, table, grouped_metadata)
      # Heatmap (middle right)
      plot_object <- reactive(mb_plot_object_shiny(display_table()$table, input, templates = templates)) %>%
        bindCache(controller()$table_plot_hash, controller()$enrichment_term, isTRUE(input$collapsed_clusters), selected_groups()) %>%
        bindEvent(display_table(), controller()$table_plot_hash, ignoreInit = FALSE, ignoreNULL = TRUE)

      mb_top_plot <- reactive(summary_track_allsamples(
        attr(table()$table, "summary_cov"),
        template = templates$cov_panel_columns_plotly
      )) %>%
        bindCache(controller()$table_hash) %>%
        bindEvent(table(), ignoreInit = FALSE, ignoreNULL = TRUE)

      mb_mid_plot <- reactive(mb_mid_plot_shiny(plot_object(), input$plotType)) %>%
        bindEvent(plot_object(), ignoreInit = FALSE, ignoreNULL = TRUE)

      mb_bottom_plot <- reactive(get_megabrowser_annotation_plot_shiny(controller, templates = templates)) %>%
        bindCache(controller()$table_hash) %>%
        bindEvent(controller()$table_hash, ignoreInit = FALSE, ignoreNULL = TRUE)

      output$myPlotlyPlot <- renderPlotly({
        req(input$plotType == "plotly")
        megabrowser_mid_reset_plot(mb_mid_plot(), controller, display_table, ns)
      }) %>%
        bindCache(controller()$table_plot_hash, controller()$enrichment_term,
                  isTRUE(input$collapsed_clusters), selected_groups(), ns("myPlotlyPlot")) %>%
        bindEvent(plot_object(), ignoreInit = FALSE, ignoreNULL = TRUE)

      mb_mid_image <- reactive(mb_mid_image_shiny(plot_object(), session, ns)) %>%
        bindEvent(plot_object(), ignoreInit = FALSE, ignoreNULL = TRUE)

      output$mb_subplot_static <- renderPlotly({
        req(input$plotType == "ggplot2")
        mb_static_subplot_shiny(mb_top_plot(), mb_mid_image(), mb_bottom_plot())
      }) %>%
        bindEvent(plot_object(), ignoreInit = FALSE, ignoreNULL = TRUE)

      output$mb_top_summary <- renderPlotly({
        megabrowser_peer_reset_plot(
          mb_top_plot(), controller, table,
          c(ns("myPlotlyPlot"), ns("mb_bottom_gene"))
        )
      }) %>%
        bindCache(controller()$table_hash) %>%
        bindEvent(mb_top_plot(), ignoreInit = FALSE, ignoreNULL = TRUE)

      output$mb_bottom_gene <- renderPlotly({
        megabrowser_peer_reset_plot(
          mb_bottom_plot(), controller, table,
          c(ns("myPlotlyPlot"), ns("mb_top_summary"))
        )
      }) %>%
        bindCache(controller()$table_hash) %>%
        bindEvent(mb_bottom_plot(), ignoreInit = FALSE, ignoreNULL = TRUE)

      rendered_plot_type <- reactiveVal(NULL)
      observeEvent(plot_object(), {
        type <- controller()$plotType
        if (!identical(isolate(rendered_plot_type()), type)) rendered_plot_type(type)
      })
      observeEvent(rendered_plot_type(), {
        shinyjs::toggleState("reset_view", condition = rendered_plot_type() == "plotly")
        shinyjs::toggleState("heatmap_only", condition = rendered_plot_type() == "plotly")
        shinyjs::toggleClass("workspace", "mega-static-view", condition = rendered_plot_type() != "plotly")
      })
      output$c <- renderUI({
        req(rendered_plot_type())
        renderMegabrowser(rendered_plot_type(), ns, height = "var(--mega-height, 700px)")
      })

      # Additional plots and tables
      observeEvent(table(), {
        choices <- c("Grouping" = "grouping", stats::setNames(names(metadata), names(metadata)))
        selected <- isolate(input$enrichment_metadata)
        if (is.null(selected) || !selected %in% choices) selected <- "grouping"
        updateSelectizeInput(session, "enrichment_metadata", choices = choices, selected = selected, server = TRUE)
      })
      enrichment_field <- reactive({
        value <- input$enrichment_metadata
        if (isTruthy(value)) value else "grouping"
      })
      meta_and_clusters <- reactive(megabrowser_metadata_enrichment(
        grouped_metadata(), metadata, enrichment_field())) %>%
        bindCache(controller()$table_hash, controller()$enrichment_term, enrichment_field())

      output$d <- renderPlotly({
        allsamples_sidebar_plotly(display_table()$meta, templates = templates)
      }) %>%
        bindCache(controller()$table_hash, controller()$enrichment_term, isTRUE(input$collapsed_clusters), selected_groups()) %>%
        bindEvent(display_table(),
                  ignoreInit = FALSE,
                  ignoreNULL = TRUE)
      output$e <- renderPlotly(allsamples_enrich_bar_plotly(meta_and_clusters()$enrich_dt)) %>%
        bindCache(controller()$table_hash, controller()$enrichment_term, enrichment_field()) %>%
        bindEvent(meta_and_clusters(),
                  ignoreInit = FALSE,
                  ignoreNULL = TRUE)

      output$stats <- renderDT(allsamples_meta_stats_shiny(meta_and_clusters()$enrich_dt)) %>%
        bindEvent(meta_and_clusters(),
                  ignoreInit = FALSE,
                  ignoreNULL = TRUE)

      module_additional_megabrowser(input, output, session)
      return(rv)
    }
  )
}
