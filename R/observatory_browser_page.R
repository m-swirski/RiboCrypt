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

#' Normalize a run vector used by observatory browser selections.
#' @noRd
observatory_clean_run_vector <- function(runs) {
  runs <- unique(as.character(runs %||% character()))
  runs[!is.na(runs) & nzchar(runs)]
}

#' Selection ids known to the observatory browser.
#' @noRd
observatory_browser_selection_ids <- function(selections, selection_index = NULL) {
  as.character(selection_index %||% names(selections) %||% character())
}

#' Active selection id constrained to the known selection ids.
#' @noRd
observatory_browser_active_selection_id <- function(active_selection_id, selection_ids) {
  selection_id <- as.character(active_selection_id %||% selection_ids[[1]])
  if (selection_id %in% selection_ids) selection_id else selection_ids[[1]]
}

#' TRUE when an empty selection should expand to all available runs.
#' @noRd
observatory_browser_should_expand_selection <- function(id, active_id, ids, label) {
  identical(label, "All merged") || (id == active_id && length(ids) == 1L)
}

#' Ensure every visible selection id has a slot in the selection list.
#' @noRd
observatory_browser_ensure_selection_slots <- function(selections, selection_ids) {
  for (id in selection_ids) {
    if (is.null(selections[[id]])) selections[[id]] <- NULL
  }
  selections
}

#' Expand empty All merged selections to the full run list.
#' @noRd
observatory_browser_expand_empty_selections <- function(selections, selection_ids,
                                                       active_id, labels, all_runs) {
  for (id in selection_ids) {
    selection <- observatory_clean_run_vector(selections[[id]])
    label <- labels[[id]] %||% ""
    if (length(selection) == 0 && length(all_runs) > 0 &&
        observatory_browser_should_expand_selection(id, active_id, selection_ids, label)) {
      selections[[id]] <- all_runs
    }
  }
  selections
}

#' Create the default observatory browser selection when none exists.
#' @noRd
observatory_browser_default_selection <- function(selections, active_selection_id, all_runs) {
  if (!length(all_runs)) return(selections)
  selection_id <- as.character(active_selection_id %||% "1")
  selections[[selection_id]] <- all_runs
  selections
}

observatory_resolve_library_selections <- function(
  library_selections,
  all_library_runs,
  active_selection_id = NULL,
  library_selection_labels = NULL,
  selection_index = NULL
) {
  selections <- library_selections %||% list()
  all_runs <- observatory_clean_run_vector(all_library_runs)

  selection_ids <- observatory_browser_selection_ids(selections, selection_index)
  if (length(selection_ids) == 0) {
    return(observatory_browser_default_selection(selections, active_selection_id, all_runs))
  }

  selections <- observatory_browser_ensure_selection_slots(selections, selection_ids)
  selection_id <- observatory_browser_active_selection_id(active_selection_id, selection_ids)
  observatory_browser_expand_empty_selections(
    selections, selection_ids, selection_id,
    library_selection_labels %||% list(), all_runs
  )
}

#' Static user info used by the observatory browser cache.
#' @noRd
observatory_browser_cache_user_info <- function() {
  list(userAgent = NULL, is_cellphone = FALSE, width = NULL, height = NULL)
}

#' Resolved library selections reactive for observatory browser plots.
#' @noRd
observatory_browser_library_selections <- function(library_selections,
                                                   all_library_runs,
                                                   active_selection_id,
                                                   library_selection_labels,
                                                   selection_index) {
  shiny::reactive({
    observatory_resolve_library_selections(
      library_selections = library_selections(),
      all_library_runs = all_library_runs(),
      active_selection_id = active_selection_id(),
      library_selection_labels = if (is.null(library_selection_labels)) NULL else library_selection_labels(),
      selection_index = if (is.null(selection_index)) NULL else selection_index()
    )
  })
}

#' Build observatory browser main-plot controller data.
#' @noRd
observatory_browser_main_controller <- function(input, tx, cds, df,
                                                library_selections,
                                                library_selection_labels) {
  click_plot_browser_main_controller(
    input = input, tx = tx, cds = cds, libs = NULL, df = df,
    user_info = function() list(is_cellphone = FALSE, width = NULL),
    library_selections = library_selections,
    library_selection_labels = library_selection_labels
  )
}

#' Plotly browser plot for observatory library selections.
#' @noRd
observatory_browser_plot <- function(main_plot_controls, bottom_panel, session,
                                     templates = NULL) {
  browser_track_panel_shiny(
    main_plot_controls, bottom_panel, session,
    ylabels = names(main_plot_controls()$library_selections),
    profiles = main_plot_controls()$profiles,
    use_fst = TRUE,
    selected_libraries = main_plot_controls()$library_selections,
    templates = templates
  )
}

#' URL context passed to shared browser auxiliary outputs.
#' @noRd
observatory_browser_url_context <- function(selected_experiment, color_by,
                                           selection_index, library_selections,
                                           library_selection_labels,
                                           active_selection_id,
                                           observatory_url_state,
                                           kickoff) {
  list(
    selected_experiment = selected_experiment,
    color_by = color_by,
    selection_index = selection_index,
    library_selections = library_selections,
    library_selection_labels = library_selection_labels,
    active_selection_id = active_selection_id,
    observatory_url_state = observatory_url_state,
    kickoff = kickoff
  )
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
    cache_user_info <- observatory_browser_cache_user_info()
    resolved_library_selections <- observatory_browser_library_selections(
      library_selections, all_library_runs, active_selection_id,
      library_selection_labels, selection_index
    )

    main_plot_controls <- shiny::reactive({
      observatory_browser_main_controller(
        input, tx, cds, df, resolved_library_selections,
        library_selection_labels
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
      observatory_browser_plot(
        main_plot_controls, bottom_panel(), session,
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
                              observatory = observatory_browser_url_context(
                                selected_experiment, color_by, selection_index,
                                library_selections, library_selection_labels,
                                active_selection_id, observatory_url_state,
                                kickoff
                              )
    )
  })
}
