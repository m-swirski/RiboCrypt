metadata_ui <- function(id, all_exp, all_exp_meta, label = "metadata") {
  ns <- NS(id)
  genomes <- unique(all_exp$organism)
  experiments <- all_exp$name
  all_merged <- all_exp[grep("all_merged", name),][libtypes == "RFP",][grep("Escherichia_coli", name, invert = TRUE),]
  navbarMenu(
    title = "Metadata", icon = icon("layer-group"),
    sample_info_ui("sample_info"),
    study_info_ui("study_info"),
    sra_search_ui("sra_search"),
    predicted_translons_ui("predicted_translons", all_merged),
    umap_ui("umap", all_exp_meta)
  )
}

metadata_server <- function(id, all_experiments, metadata, all_exp_meta,
                            browser_options) {
  if (!is.null(metadata)) {
    sample_info_server("sample_info", metadata)
  } else print("No metadata given, ignoring Sample_info server.")
  study_info_server("study_info", all_experiments)
  sra_search_server("sra_search")
  predicted_translons_server("predicted_translons", all_experiments, browser_options)
  umap_server("umap", all_exp_meta, browser_options)
}
