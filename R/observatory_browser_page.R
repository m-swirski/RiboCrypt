observatory_browser_ui <- function(id, browser_options) {
  browser_ui_shared(
    id = id,
    browser_options = browser_options,
    title = "Browse",
    icon_name = "chart-line",
    include_global_init = FALSE,
    include_browser_sources = FALSE,
    include_library_controls = FALSE,
    include_download_button = FALSE,
    include_aux_outputs = FALSE,
    main_plot_output_id = "browser_plot",
    main_plot_height = "700px"
  )
}

# TODO
# Add settings for kmer, leader and trailer extensions and aggregation method
observatory_browser_server <- function(
  id,
  df,
  library_selections,
  library_selection_labels,
  gene_name_list, tx, cds, experiments, org, gg_theme,
  rv, browser_options,
  selection_index = shiny::reactive(NULL),
  active_selection_id = shiny::reactive(NULL),
  selected_experiment = shiny::reactive(NULL),
  color_by = shiny::reactive(NULL),
  observatory_url_state = shiny::reactiveVal(NULL)
) {
  shiny::moduleServer(id, function(input, output, session) {
    # -- Gene / transcript input wiring ------------------------------------
    allsamples_observer_controller(input, output, session)
    kickoff <- shiny::reactiveVal(FALSE)
    module_browser_shared_ui(input, output, session, function() {
      state <- observatory_state_from_inputs(
        selected_experiment = selected_experiment(),
        color_by = color_by(),
        view = "browser",
        browser = list(
          gene = input$gene,
          tx = input$tx,
          frames_type = input$frames_type,
          kmer = input$kmer,
          extendLeaders = input$extendLeaders,
          extendTrailers = input$extendTrailers,
          viewMode = input$viewMode,
          other_tx = input$other_tx,
          collapsed_introns = input$collapsed_introns,
          collapsed_introns_width = input$collapsed_introns_width,
          genomic_region = input$genomic_region,
          zoom_range = input$zoom_range,
          customSequence = input$customSequence,
          go = TRUE
        ),
        selections = list(
          index = selection_index(),
          plot_selections = library_selections(),
          data_table_selections = library_selections(),
          labels = library_selection_labels(),
          active_selection_id = active_selection_id()
        )
      )
      observatory_clipboard_url_button(state, session)
    })

    shiny::observe({
      st <- observatory_url_state()
      if (is.null(st)) return()
      br <- st$browser
      if (is.null(br)) return()

      if (!is.null(br$gene)) shiny::updateSelectizeInput(session, "gene", selected = br$gene)
      if (!is.null(br$tx)) shiny::updateSelectizeInput(session, "tx", selected = br$tx)
      if (!is.null(br$frames_type)) shiny::updateSelectizeInput(session, "frames_type", selected = br$frames_type)
      if (!is.null(br$kmer) && !is.na(as.numeric(br$kmer))) shiny::updateSliderInput(session, "kmer", value = as.numeric(br$kmer))
      if (!is.null(br$extendLeaders) && !is.na(as.numeric(br$extendLeaders))) {
        shiny::updateNumericInput(session, "extendLeaders", value = as.numeric(br$extendLeaders))
      }
      if (!is.null(br$extendTrailers) && !is.na(as.numeric(br$extendTrailers))) {
        shiny::updateNumericInput(session, "extendTrailers", value = as.numeric(br$extendTrailers))
      }
      if (!is.null(br$viewMode)) shinyWidgets::updatePrettySwitch(session, "viewMode", value = isTRUE(br$viewMode))
      if (!is.null(br$other_tx)) shinyWidgets::updatePrettySwitch(session, "other_tx", value = isTRUE(br$other_tx))
      if (!is.null(br$collapsed_introns)) {
        shinyWidgets::updatePrettySwitch(session, "collapsed_introns", value = isTRUE(br$collapsed_introns))
      }
      if (!is.null(br$collapsed_introns_width) && !is.na(as.numeric(br$collapsed_introns_width))) {
        shiny::updateNumericInput(session, "collapsed_introns_width", value = as.numeric(br$collapsed_introns_width))
      }
      if (!is.null(br$genomic_region)) shiny::updateTextInput(session, "genomic_region", value = br$genomic_region)
      if (!is.null(br$zoom_range)) shiny::updateTextInput(session, "zoom_range", value = br$zoom_range)
      if (!is.null(br$customSequence)) shiny::updateTextInput(session, "customSequence", value = br$customSequence)
    }) |> shiny::bindEvent(observatory_url_state(), ignoreInit = FALSE, once = TRUE)

    shiny::observe({
      st <- observatory_url_state()
      if (is.null(st) || !identical(st$view, "browser")) return()
      if (!isTRUE(st$browser$go)) return()
      req(nzchar(input$gene), nzchar(input$tx))
      sel <- library_selections()
      if (is.null(sel) || !any(lengths(sel) > 0)) return()
      kickoff(TRUE)
    }) |> shiny::bindEvent(
      input$gene, input$tx,
      library_selections(),
      observatory_url_state(),
      ignoreInit = TRUE
    )


    main_plot_controls <- shiny::reactive({
      click_plot_browser_main_controller(
        input = input,
        tx = tx,
        cds = cds,
        libs = NULL,
        df = df,
        gg_theme = gg_theme,
        user_info = function() list(is_cellphone = FALSE, width = NULL),
        library_selections = library_selections,
        library_selection_labels = library_selection_labels
      )
    }) |> shiny::bindEvent(list(input$go, kickoff()), ignoreNULL = TRUE, ignoreInit = TRUE)

    # -- Plot rendering ----------------------------------------------------

    bottom_panel <- shiny::reactive({
      bottom_panel_shiny(main_plot_controls)
    }) |>
      shiny::bindCache(main_plot_controls()$hash_bottom) |>
      shiny::bindEvent(main_plot_controls(), ignoreNULL = TRUE)

    browser_plot <- shiny::reactive({
      browser_track_panel_shiny(
        main_plot_controls, bottom_panel(), session,
        ylabels = names(main_plot_controls()$library_selections),
        profiles = main_plot_controls()$profiles,
        use_fst = TRUE,
        selected_libraries = main_plot_controls()$library_selections
      )
    }) |> shiny::bindEvent(bottom_panel(), ignoreNULL = TRUE)

    output$browser_plot <- plotly::renderPlotly({
      browser_plot()
    }) |> shiny::bindEvent(browser_plot(), ignoreNULL = TRUE)
  })
}
