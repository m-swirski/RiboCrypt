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

observatory_resolve_library_selections <- function(
  library_selections,
  all_library_runs,
  active_selection_id = NULL,
  library_selection_labels = NULL,
  selection_index = NULL
) {
  selections <- library_selections %||% list()
  all_runs <- unique(as.character(all_library_runs %||% character()))
  all_runs <- all_runs[!is.na(all_runs) & nzchar(all_runs)]

  selection_ids <- as.character(selection_index %||% names(selections) %||% character())
  if (length(selection_ids) == 0) {
    if (length(all_runs) == 0) return(selections)
    selection_id <- as.character(active_selection_id %||% "1")
    selections[[selection_id]] <- all_runs
    return(selections)
  }

  for (id in selection_ids) {
    if (is.null(selections[[id]])) selections[[id]] <- NULL
  }

  selection_id <- as.character(active_selection_id %||% selection_ids[[1]])
  if (!(selection_id %in% selection_ids)) selection_id <- selection_ids[[1]]

  for (id in selection_ids) {
    selection <- as.character(selections[[id]] %||% character())
    selection <- unique(selection[!is.na(selection) & nzchar(selection)])
    label <- library_selection_labels[[id]] %||% ""
    should_expand_to_all <- identical(label, "All merged") || (id == selection_id && length(selection_ids) == 1L)
    if (length(selection) == 0 && length(all_runs) > 0 && should_expand_to_all) {
      selections[[id]] <- all_runs
    }
  }

  selections
}

# TODO
# Add settings for kmer, leader and trailer extensions and aggregation method
observatory_browser_server <- function(
  id,
  df,
  library_selections,
  all_library_runs,
  library_selection_labels,
  gene_name_list, tx, cds, experiments, org,
  rv, browser_options,
  templates = NULL,
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
    cache_user_info <- list(userAgent = NULL, is_cellphone = FALSE, width = NULL, height = NULL)

    main_plot_controls <- shiny::reactive({
      click_plot_browser_main_controller(
        input = input,
        tx = tx,
        cds = cds,
        libs = NULL,
        df = df,
        user_info = function() list(is_cellphone = FALSE, width = NULL),
        library_selections = shiny::reactive({
          observatory_resolve_library_selections(
            library_selections = library_selections(),
            all_library_runs = all_library_runs(),
            active_selection_id = active_selection_id(),
            library_selection_labels = if (is.null(library_selection_labels)) NULL else library_selection_labels(),
            selection_index = if (is.null(selection_index)) NULL else selection_index()
          )
        }),
        library_selection_labels = library_selection_labels
      )
    }) |>
      shiny::bindCache(
        input_to_list(input, cache_user_info),
        observatory_selection_cache_key(
          library_selections(),
          if (is.null(library_selection_labels)) NULL else library_selection_labels()
        ),
        selected_experiment()
      ) |>
      shiny::bindEvent(list(input$go, kickoff()), ignoreNULL = TRUE, ignoreInit = TRUE)

    # -- Plot rendering ----------------------------------------------------

    bottom_panel <- shiny::reactive({
      bottom_panel_shiny(main_plot_controls, templates = templates)
    }) |>
      shiny::bindCache(main_plot_controls()$hash_bottom) |>
      shiny::bindEvent(main_plot_controls(), ignoreNULL = TRUE)

    browser_plot <- shiny::reactive({
      browser_track_panel_shiny(
        main_plot_controls, bottom_panel(), session,
        ylabels = names(main_plot_controls()$library_selections),
        profiles = main_plot_controls()$profiles,
        use_fst = TRUE,
        selected_libraries = main_plot_controls()$library_selections,
        templates = templates
      )
    }) |>
      shiny::bindCache(main_plot_controls()$hash_browser) |>
      shiny::bindEvent(bottom_panel(), ignoreNULL = TRUE)

    output$browser_plot <- plotly::renderPlotly({
      browser_plot()
    }) |>
      shiny::bindCache(main_plot_controls()$hash_browser) |>
      shiny::bindEvent(browser_plot(), ignoreNULL = TRUE)

    module_additional_browser(input, output, session,
                              mode = "observatory",
                              observatory = list(
                                selected_experiment = selected_experiment,
                                color_by = color_by,
                                selection_index = selection_index,
                                library_selections = library_selections,
                                library_selection_labels = library_selection_labels,
                                active_selection_id = active_selection_id,
                                observatory_url_state = observatory_url_state,
                                kickoff = kickoff
                              )
    )
  })
}
