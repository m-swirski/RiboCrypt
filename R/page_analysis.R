analysis_ui <- function(id, all_exp, browser_options, libs, metadata, all_exp_meta,
                        label = "Analysis") {
  ns <- NS(id)
  genomes <- unique(all_exp$organism)
  experiments <- all_exp$name
  navbarMenu(
    title = "Analysis", icon = icon("layer-group"),
    heatmap_ui("Heatmap", all_exp, browser_options, libs),
    codon_ui("Codon", all_exp, browser_options, libs),
    DEG_ui("DEG", all_exp, browser_options),
    quality_ui("quality", all_exp, browser_options, libs),
    fastq_ui("FastQ_report", all_exp, browser_options, libs)
  )
}

analysis_server <- function(id, all_experiments, without_readlengths_env,
                         with_readlengths_env, df, df_with, experiments,
                         tx, cds, libs, org, gene_name_list, rv, metadata,
                         all_exp_meta, exp_init_meta, df_meta, names_init, browser_options) {
  rv <- heatmap_server("Heatmap", all_experiments, with_readlengths_env,
                df_with, experiments, tx, cds, libs, org, gene_name_list, rv)
  rv <- codon_server("Codon", all_experiments, without_readlengths_env,
                     df, experiments, tx, cds, libs, org, gene_name_list, rv)
  rv <- DEG_server("DEG", all_experiments, without_readlengths_env, df,
                   experiments, libs, org, gene_name_list, rv)
  rv <- quality_server("quality", all_experiments, with_readlengths_env,
                df_with, experiments, tx, cds, libs, org, gene_name_list, rv)
  rv <- fastq_server("FastQ_report", all_experiments, df, experiments, libs, org, rv)
  return(rv)
}
