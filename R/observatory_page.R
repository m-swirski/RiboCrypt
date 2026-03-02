observatory_ui <- function(id, meta_experiment_list, browser_options) {
  ns <- shiny::NS(id)

  shiny::tabPanel(
    "Observatory",
    shiny::tabsetPanel(
      observatory_selector_ui(ns("selector"), meta_experiment_list, browser_options),
      observatory_browser_ui(ns("browser"))
    )
  )
}

observatory_server <- function(
  id,
  meta_experiment_list, all_libraries_df, names_init_meta
) {
  shiny::moduleServer(id, function(input, output, session) {
    selections <- observatory_selector_server(
      "selector",
      meta_experiment_list,
      all_libraries_df
    )

    observatory_browser_server(
      "browser",
      selections$meta_experiment_df,
      selections$selected_libraries$data_table_selections,
      selections$selected_libraries$labels
    )
  })
}
