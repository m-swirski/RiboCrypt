analysis_ui <- function(id, all_exp, browser_options, libs, label = "Analysis") {
  ns <- NS(id)
  genomes <- unique(all_exp$organism)
  experiments <- all_exp$name
  navbarMenu(
    title = "analysis", icon = icon("layer-group"),
    heatmap_ui("heatmap", all_exp, browser_options, libs),
    codon_ui("codon", all_exp, browser_options, libs),
    DEG_ui("DEG", all_exp, browser_options),
    quality_ui("quality", all_exp, browser_options, libs),
    fastq_ui("fastq", all_exp, browser_options, libs),
    browser_allsamp_ui("browser_allsamp", all_exp, browser_options, libs)
  )
}

analysis_server <- function(id, all_experiments, without_readlengths_env,
                         with_readlengths_env, df, df_with, experiments,
                         tx, cds, libs, org, gene_name_list, rv) {
  rv <- heatmap_server("heatmap", all_experiments, with_readlengths_env,
                df_with, experiments, tx, cds, libs, org, gene_name_list, rv)
  rv <- codon_server("codon", all_experiments, without_readlengths_env,
                     df, experiments, tx, cds, libs, org, gene_name_list, rv)
  rv <- DEG_server("DEG", all_experiments, without_readlengths_env, df,
                   experiments, libs, org, rv)
  rv <- quality_server("quality", all_experiments, with_readlengths_env,
                df_with, experiments, tx, cds, libs, org, gene_name_list, rv)
  rv <- fastq_server("fastq", all_experiments, df, experiments, libs, org, rv)
  browser_allsamp_server("browser_allsamp", all_experiments, without_readlengths_env, df,
                         experiments, tx, cds, libs, org, gene_name_list, rv)
  return(rv)
}
