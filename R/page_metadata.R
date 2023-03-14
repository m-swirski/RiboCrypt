metadata_ui <- function(id, all_exp, label = "metadata") {
  ns <- NS(id)
  genomes <- unique(all_exp$organism)
  experiments <- all_exp$name
  navbarMenu(
    title = "metadata", icon = icon("layer-group"),
    sra_search_ui("sra_search"),
    study_info_ui("study_info")
  )
}

metadata_server <- function(id, all_experiments) {
  sra_search_server("sra_search")
  study_info_server("study_info", all_experiments)
}
