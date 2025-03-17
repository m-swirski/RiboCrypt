metadata_ui <- function(id, all_exp, metadata, label = "metadata") {
  ns <- NS(id)
  genomes <- unique(all_exp$organism)
  experiments <- all_exp$name
  all_merged <- all_exp[grep("all_merged", name),][libtypes == "RFP",][grep("Escherichia_coli", name, invert = TRUE),]
  navbarMenu(
    title = "Metadata", icon = icon("layer-group"),
    run_info_ui("run_info"),
    study_info_ui("study_info"),
    sra_search_ui("sra_search"),
    predicted_translons_ui("predicted_translons", all_merged)
  )
}

metadata_server <- function(id, all_experiments, metadata) {
  run_info_server("run_info", metadata)
  study_info_server("study_info", all_experiments)
  sra_search_server("sra_search")
  predicted_translons_server("predicted_translons")
}
