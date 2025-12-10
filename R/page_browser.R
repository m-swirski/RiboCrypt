browser_ui <- function(id, all_exp, browser_options, gene_names_init,
                       libs, label = "Browser") {

  ns <- NS(id)
  genomes <- unique(all_exp$organism)
  experiments <- all_exp$name
  all_isoforms <- subset(gene_names_init, label == browser_options["default_gene"])
  init_libs <- unlist(strsplit(browser_options["default_libs"], "\\|"))
  viewMode <- browser_options["default_view_mode"] == "genomic"
  introns_width <- as.numeric(browser_options["collapsed_introns_width"])
  full_annotation <- as.logical(browser_options["full_annotation"])
  translons <- as.logical(browser_options["translons"])
  normalizations <- normalizations("metabrowser")
  panel_hidden_or_not_class <- ifelse(browser_options["hide_settings"] == "TRUE",
                                      "floating_settings_panel hidden",
                                      "floating_settings_panel")
  tabPanel(
    title = "browser", icon = icon("chart-line"),
    shinyjs::useShinyjs(),
    rclipboardSetup(),
    tags$head(includeHTML(system.file("google_analytics_html", "google_analytics.html", package = "RiboCrypt"))),
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
              tabPanel("Browser",
                       fluidRow(
                         column(6, organism_input_select(c("ALL", genomes), ns)),
                         column(6, experiment_input_select(experiments, ns, browser_options))
                       ),
                       fluidRow(
                         column(11, library_input_select(ns, TRUE, libs, init_libs)),
                         column(1, actionButton(ns("select_all_btn"), "", icon = icon("check"),
                                                class = "btn btn-sm btn-primary", title = "Select all"))
                       ),
                       fluidRow(
                         column(6, frame_type_select(ns, selected = browser_options["default_frame_type"])),
                         column(6, sliderInput(ns("kmer"), "K-mer length", min = 1, max = 20,
                                               value = as.numeric(browser_options["default_kmer"])))
                       ),
                       fluidRow(
                         column(6, normalization_input_select(ns, choices = normalizations,
                                                              help_link = "mbrowser"))
                       )
              ),
              tabPanel("Settings",
                                fluidRow(
                                  numericInput(ns("extendLeaders"), "5' extension", 0),
                                  numericInput(ns("extendTrailers"), "3' extension", 0),
                                  numericInput(ns("collapsed_introns_width"), "Collapse Introns (nt flanks)",
                                               introns_width)
                                ),
                                fluidRow(textInput(ns("genomic_region"), "Genomic region", ""),
                                         textInput(ns("zoom_range"), "Zoom interval", ""),
                                         textInput(ns("customSequence"), "Custom sequences highlight", "")
                                ),
                                fluidRow(
                                  checkboxInput(ns("add_uorfs"), "uORF annotation", FALSE),
                                  checkboxInput(ns("add_translon"), "Predicted translons", translons),
                                  ),
                                fluidRow(checkboxInput(ns("log_scale"), "Log scale", FALSE),
                                         checkboxInput(ns("log_scale_protein"), "Log scale Protein", FALSE)
                                ),
                                fluidRow(
                                  column(4, checkboxInput(ns("expression_plot"), "Gene expression plot", FALSE)),
                                  column(4, checkboxInput(ns("useCustomRegions"), "Protein structures", TRUE)),
                                ),
                                fluidRow(column(4, checkboxInput(ns("phyloP"), "Conservation (phyloP)", FALSE)),
                                         column(4, checkboxInput(ns("mapability"), "Mapability (28mers)", FALSE))),
                                fluidRow(
                                  column(6, checkboxInput(ns("withFrames"), "Split color Frames", TRUE)),
                                  column(6, frame_subsetter_select(ns))
                                ),
                                fluidRow(
                                  column(6, checkboxInput(ns("summary_track"), "Summary top track", FALSE)),
                                  column(6, frame_type_select(ns, "summary_track_type", "Summary display type"))
                                ),
                                fluidRow(
                                  column(4, downloadButton(ns("download_plot_html"), "Download HTML",
                                               style = "width: 100%; font-size: 14px; font-weight: bold; background-color: #007bff; color: white; border-color: white !important;")),
                                  column(4, export_format_of_plot(ns)),
                                  column(4, uiOutput(ns("clip"))))
              )
            )
      ))),
      column(2, gene_input_select(ns, FALSE, browser_options)),
      column(2, tx_input_select(ns, FALSE, all_isoforms, browser_options["default_isoform"])),
      column(1, NULL, plot_button(ns("go"))),
      column(1, prettySwitch(ns("viewMode"), "Genomic View", value = viewMode,
                             status = "success", fill = TRUE, bigger = TRUE),
                prettySwitch(ns("other_tx"), "Full annotation", value = full_annotation,
                              status = "success", fill = TRUE, bigger = TRUE),
                prettySwitch(ns("collapsed_introns"), "Collapse introns", value = FALSE,
                             status = "success", fill = TRUE, bigger = TRUE))
    ),
    tags$hr(),
    # ---- Full Width Main Panel ----
    fluidRow(
      column(12,
             jqui_resizable(plotlyOutput(ns("c"), height = "500px")) %>% shinycssloaders::withSpinner(color="#0dc5c1"),
             plotlyOutput(ns("e"), height = "50px"),
             uiOutput(ns("proteinStruct")),
             plotlyOutput(ns("d")) %>% shinycssloaders::withSpinner(color="#0dc5c1")
      )
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
      mainPlotControls <- reactive({
        click_plot_browser_main_controller(
          tx, 
          cds, 
          libs, 
          df,
          selectedGene = input$gene, 
          selectedTx = input$tx, 
          otherTx = input$other_tx, 
          addUorfs = input$add_uorfs,
          addTranslons = input$add_translon,
          collapsedIntrons = input$collapsed_introns,
          collapsedIntronsWidth = input$collapsed_introns_width,
          genomicRegion = input$genomic_region,
          extendLeaders = input$extendLeaders,
          extendTrailers = input$extendTrailers,
          zoomRange = input$zoom_range,
          viewMode = input$viewMode,
          selectedLibraries = input$library,
          withFrames = input$withFrames,
          exportFormat = input$plot_export_format,
          summaryTrack = input$summary_track,
          summaryTrackType = input$summary_track_type,
          kmer = input$kmer,
          framesType = input$frames_type,
          framesSubset = input$frames_subset,
          normalization = input$normalization,
          customSequence = input$customSequence,
          logScale = input$log_scale,
          phyloP = input$phyloP,
          mapability = input$mapability,
          expressionPlot = input$expression_plot
          )
        }) %>% bindEvent(input$go, ignoreInit = check_plot_on_start(browser_options), ignoreNULL = FALSE)

      bottom_panel <- reactive(bottom_panel_shiny(mainPlotControls))  %>%
        bindCache(mainPlotControls()$hash_bottom) %>%
        bindEvent(mainPlotControls(), ignoreInit = FALSE, ignoreNULL = TRUE)

      browser_plot <- reactive(browser_track_panel_shiny(mainPlotControls, bottom_panel(), session)) %>%
        bindCache(mainPlotControls()$hash_browser) %>%
        bindEvent(bottom_panel(), ignoreInit = FALSE, ignoreNULL = TRUE)

      output$c <- renderPlotly(browser_plot()) %>%
        bindCache(mainPlotControls()$hash_browser) %>%
        bindEvent(browser_plot(), ignoreInit = FALSE, ignoreNULL = TRUE)

      # Additional outputs
      module_additional_browser(input, output, session)

      return(rv)
    }
  )
}
