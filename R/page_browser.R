browser_ui_shared <- function(id, browser_options, gene_names_init = NULL,
                              all_exp = NULL, libs = NULL,
                              label = "browser",
                              title = "browser",
                              icon_name = "chart-line",
                              include_global_init = TRUE,
                              include_browser_sources = TRUE,
                              include_library_controls = TRUE,
                              include_download_button = TRUE,
                              include_aux_outputs = TRUE,
                              main_plot_output_id = "c",
                              main_plot_height = "550px") {
  browser_option <- function(name, default = "") {
    value <- browser_options[name]
    if (is.null(value) || length(value) == 0 || is.na(value)) default else value[[1]]
  }

  ns <- NS(id)
  genomes <- if (include_browser_sources) unique(all_exp$organism) else character()
  experiments <- if (include_browser_sources) all_exp$name else character()
  all_isoforms <- if (is.null(gene_names_init)) {
    data.frame(label = character(), stringsAsFactors = FALSE)
  } else {
    subset(gene_names_init, label == browser_options["default_gene"])
  }
  init_libs <- if (include_library_controls) {
    unlist(strsplit(browser_option("default_libs"), "\\|"))
  } else {
    character()
  }
  viewMode <- identical(browser_option("default_view_mode", "tx"), "genomic")
  introns_width <- suppressWarnings(as.numeric(browser_option("collapsed_introns_width", "100")))
  if (is.na(introns_width)) introns_width <- 100
  full_annotation <- isTRUE(as.logical(browser_option("full_annotation", "FALSE")))
  translons <- isTRUE(as.logical(browser_option("translons", "FALSE")))
  translons_transcode <- isTRUE(as.logical(browser_option("translons_transcode", "FALSE")))
  panel_hidden_or_not_class <- ifelse(browser_option("hide_settings", "TRUE") == "TRUE",
                                      "floating_settings_panel hidden",
                                      "floating_settings_panel")
  browser_tab_children <- list(
    if (include_global_init) shinyjs::useShinyjs(),
    if (include_global_init) rclipboardSetup(),
    if (include_global_init) tags$head(includeHTML(system.file("google_analytics_html", "google_analytics.html", package = "RiboCrypt"))),
    if (include_global_init) tags$script(HTML(sprintf("
      $(document).on('shiny:connected', function() {
        Shiny.setInputValue('%s', navigator.userAgent, {priority: 'event'});
      });
    ", ns("js_user_agent")))),
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
                       if (include_browser_sources) {
                         fluidRow(
                           column(6, organism_input_select(c("ALL", genomes), ns)),
                           column(6, experiment_input_select(experiments, ns, browser_options))
                         )
                       },
                       if (include_library_controls) {
                         fluidRow(
                           column(11, library_input_select(ns, TRUE, libs, init_libs)),
                           column(1, actionButton(ns("select_all_btn"), "", icon = icon("check"),
                                                  class = "btn btn-sm btn-primary", title = "Select all"))
                         )
                       },
                       if (include_library_controls) {
                         fluidRow(prettySwitch(ns("unique_align"), "Unique alignments", value = FALSE,
                                      status = "success", fill = TRUE, bigger = TRUE))
                       },
                       fluidRow(
                         column(6, frame_type_select(ns, selected = browser_option("default_frame_type", "columns"))),
                         column(6, sliderInput(ns("kmer"), "K-mer length", min = 1, max = 20,
                                               value = as.numeric(browser_option("default_kmer", "1"))))
                       ),
                       tags$hr(style = "padding-top: 50px; padding-bottom: 50px;")
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
                                fluidRow(checkboxInput(ns("add_uorfs"), tagList("uORF annotation", tags$br(), "(all candidates)"), FALSE)),
                                fluidRow(
                                  column(4, checkboxInput(ns("add_translon"), "Predicted translons (Our all-merged: T)", translons)),
                                  column(4, checkboxInput(ns("add_translons_transcode"), "Predicted translons (TransCode: TC)",
                                                          translons_transcode))
                                  ),
                                fluidRow(
                                  column(4, checkboxInput(ns("log_scale"), "Log scale", FALSE)),
                                  column(4, checkboxInput(ns("log_scale_protein"), "Log scale Protein", FALSE))
                                ),
                                fluidRow(
                                  column(4, checkboxInput(ns("expression_plot"), "Gene expression plot", FALSE)),
                                  column(4, checkboxInput(ns("useCustomRegions"), "Protein structures", TRUE)),
                                ),
                                fluidRow(column(4, checkboxInput(ns("phyloP"), "Conservation (phyloP)", FALSE)),
                                         column(4, checkboxInput(ns("mapability"), "Mapability (28mers)", FALSE))),
                                fluidRow(
                                  column(4, checkboxInput(ns("withFrames"), "Split color Frames", TRUE)),
                                  column(4, selectizeInput(
                                    inputId = ns("colors"),
                                    label = "Frame Color theme",
                                    choices = c("R", "Color_blind")
                                  )),
                                  column(4, frame_subsetter_select(ns))
                                ),
                                fluidRow(
                                  column(4, checkboxInput(ns("summary_track"), "Summary top track", FALSE)),
                                  column(4, NULL),
                                  column(4, frame_type_select(ns, "summary_track_type", "Summary display type"))
                                ),
                                fluidRow(
                                  column(4, if (include_download_button) {
                                    downloadButton(ns("download_plot_html"), "Download HTML",
                                                   style = "width: 100%; font-size: 14px; font-weight: bold; background-color: #007bff; color: white; border-color: white !important;")
                                  }),
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
             jqui_resizable(plotlyOutput(ns(main_plot_output_id), height = main_plot_height)) %>% shinycssloaders::withSpinner(color="#0dc5c1"),
             if (include_aux_outputs) plotlyOutput(ns("e"), height = "50px"),
             if (include_aux_outputs) uiOutput(ns("proteinStruct")),
             if (include_aux_outputs) plotlyOutput(ns("d")) %>% shinycssloaders::withSpinner(color="#0dc5c1")
      )
    )
  )
  do.call(tabPanel, c(list(title = title, icon = icon(icon_name)), browser_tab_children))
}

browser_ui <- function(id, all_exp, browser_options, gene_names_init,
                       libs, label = "Browser") {
  browser_ui_shared(
    id = id,
    browser_options = browser_options,
    gene_names_init = gene_names_init,
    all_exp = all_exp,
    libs = libs,
    label = label,
    title = "browser",
    icon_name = "chart-line",
    include_global_init = TRUE,
    include_browser_sources = TRUE,
    include_library_controls = TRUE,
    include_download_button = TRUE,
    include_aux_outputs = TRUE,
    main_plot_output_id = "c",
    main_plot_height = "550px"
  )
}

browser_server <- function(id, all_experiments, env, df, experiments,
                           tx, cds, libs, org, gene_name_list, rv,
                           browser_options,
                           templates = NULL) {
  moduleServer(
    id,
    function(input, output, session, all_exp = all_experiments) {
      study_and_gene_observers(input, output, session)

      # Main plot controller, this code is only run if 'plot' is pressed
      mainPlotControls <- reactive(click_plot_browser_main_controller(input, tx, cds, libs, df, user_info)) %>%
        bindCache(input_to_list(input, user_info())) %>%
        bindEvent(list(input$go, kickoff()), ignoreInit = TRUE, ignoreNULL = FALSE)

      bottom_panel <- reactive(bottom_panel_shiny(mainPlotControls, templates = templates))  %>%
        bindCache(mainPlotControls()$hash_bottom) %>%
        bindEvent(mainPlotControls()$hash_bottom, ignoreInit = FALSE, ignoreNULL = TRUE)

      browser_plot <- reactive(browser_track_panel_shiny(
        mainPlotControls, bottom_panel(), session, templates = templates
      )) %>%
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

#' When input is ready, start ploting if specified.
#' @noRd
go_when_input_is_ready <- function(input, browser_options, fired, kickoff, libs) {
  if (fired()) return()
  if (!isTRUE(as.logical(browser_options[["plot_on_start"]]))) {
    fired(TRUE)
    return()
  }
  if (!nzchar(input$gene) || !nzchar(input$tx)) return()
  if (!identical(input$gene, browser_options[["default_gene"]])) return()
  if (!identical(input$tx,   browser_options[["default_isoform"]])) return()
  libs_wanted <- libraries_string_split(browser_options[["default_libs"]], isolate(libs()))
  if (!identical(input$library, libs_wanted)) {
    message("Libraries wanted not matching yet!")
    print(libs_wanted)
    print(isolate(input$library))
    return()
  }
  fired(TRUE)
  kickoff(TRUE)
}

#' Convert input to list, and remove irrelevant items and add user_info
#' @noRd
input_to_list <- function(input, user_info = NULL) {

  a <- reactiveValuesToList(input, TRUE)
  to_ignore <- c("go", "toggle_settings", "select_all_btn")
  resize_input_pattern <- "(__shinyjquiBookmarkState__resizable|_is_resizing|_size)$"

  a <- a[-which(names(a) %in% to_ignore | grepl(resize_input_pattern, names(a)))]
  a <- c(a, user.info = user_info[-1])
  return(a)
}
