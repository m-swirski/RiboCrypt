metadata_ui <- function(id, all_exp, all_exp_meta, label = "metadata") {
  ns <- NS(id)
  genomes <- unique(all_exp$organism)
  experiments <- all_exp$name
  navbarMenu(
    title = "Metadata", icon = icon("layer-group"),
    sample_info_ui("sample_info"),
    study_info_ui("study_info"),
    sra_search_ui("sra_search"),
    umap_ui("umap", all_exp_meta)
  )
}

metadata_server <- function(id, all_experiments, metadata, all_exp_meta,
                            browser_options) {
  if (!is.null(metadata)) {
    sample_info_server("sample_info", metadata, browser_options["search_on_init"])
  } else print("No metadata given, ignoring Sample_info server.")
  study_info_server("study_info", all_experiments)
  sra_search_server("sra_search")
  umap_server("umap", all_exp_meta, browser_options)
}
