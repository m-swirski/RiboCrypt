observatory_ui <- function(id, meta_experiment_list, browser_options) {
  ns <- shiny::NS(id)

  shiny::tabPanel(
    "Observatory",
    shiny::tabsetPanel(
      observatory_selector_ui(ns("selector"), meta_experiment_list, browser_options),
      observatory_browser_ui(ns("browser_obs"))
    )
  )
}

observatory_server <- function(
  id,
  all_exp, df, experiments,
  gene_name_list, org,
  metadata, browser_options, rv
) {
  shiny::moduleServer(id, function(input, output, session) {
    selections <- observatory_selector_server(
      "selector",
      all_exp, df, metadata, experiments, org, rv, browser_options
    )

    observatory_browser_server(
      "browser_obs",
      selections$meta_experiment_df,
      selections$selected_libraries$data_table_selections,
      selections$selected_libraries$labels, gene_name_list,
      experiments, org, rv, browser_options
    )
    return(rv)
  })
}
