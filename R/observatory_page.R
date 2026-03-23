observatory_ui <- function(id, meta_experiment_list, browser_options) {
  ns <- shiny::NS(id)

  shiny::tabPanel(
    "Observatory",
    shiny::tabsetPanel(
      id = ns("observatory_view"),
      observatory_selector_ui(ns("selector"), meta_experiment_list, browser_options),
      observatory_browser_ui(ns("browser_obs"), browser_options)
    )
  )
}

observatory_server <- function(
  id,
  all_exp, df, experiments,
  gene_name_list, tx, cds, org,
  metadata, browser_options, rv,
  templates = NULL
) {
  shiny::moduleServer(id, function(input, output, session) {
    observatory_url_state <- shiny::reactiveVal(
      parse_observatory_url(
        query = shiny::isolate(shiny::getQueryString()),
        hash = shiny::isolate(session$clientData$url_hash)
      )
    )

    shiny::observe({
      st <- observatory_url_state()
      if (is.null(st$view)) return()
      selected <- if (identical(st$view, "browser")) "Browse" else "Select libraries"
      shiny::updateTabsetPanel(session, "observatory_view", selected = selected)
    }) |> shiny::bindEvent(observatory_url_state(), ignoreInit = FALSE, once = TRUE)

    selections <- observatory_selector_server(
      "selector",
      all_exp, df, metadata, experiments, org, rv, browser_options,
      observatory_url_state = observatory_url_state
    )

    observatory_browser_server(
      "browser_obs",
      selections$meta_experiment_df,
      selections$selected_libraries$data_table_selections,
      selections$selected_libraries$labels, gene_name_list,
      tx, cds, experiments, org, rv, browser_options,
      templates = templates,
      selection_index = selections$selected_libraries$index,
      active_selection_id = selections$active_selection_id,
      selected_experiment = selections$selected_experiment,
      color_by = selections$color_by,
      observatory_url_state = observatory_url_state
    )
    return(rv)
  })
}
