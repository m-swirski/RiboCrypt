analysis_ui <- function(id, label = "Analysis", all_exp) {
  ns <- NS(id)
  genomes <- unique(all_exp$organism)
  experiments <- all_exp$name
  navbarMenu(
    title = "analysis", icon = icon("layer-group"),
    heatmap_ui("heatmap", all_exp = all_exp),
    codon_ui("codon", all_exp = all_exp),
    quality_ui("quality", all_exp = all_exp),
    fastq_ui("fastq", all_exp = all_exp)
  )
}

analysis_server <- function(id, all_experiments, without_readlengths_env,
                         with_readlengths_env) {
  heatmap_server("heatmap", all_experiments, with_readlengths_env)
  codon_server("codon", all_experiments, without_readlengths_env)
  quality_server("quality", all_experiments, with_readlengths_env)
  fastq_server("fastq", all_experiments)
}
