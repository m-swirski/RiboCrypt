observatory_ui <- function(id, meta_experiment_list) {
  ns <- shiny::NS(id)

  shiny::tabPanel(
    "Observatory",
    shiny::tabsetPanel(
      observatory_selector_ui(ns("selector"), meta_experiment_list),
      observatory_browser_ui(ns("browser"))
    )
  )
}

observatory_server <- function(
  id,
  meta_experiment_list,
  all_libraries_df
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
      selections$selected_libraries
    )
  })
}
