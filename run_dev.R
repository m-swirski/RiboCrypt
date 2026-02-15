#!/usr/bin/env Rscript

devtools::load_all(".")

all_exp <- list.experiments(
  validate = FALSE,
  BPPARAM = BiocParallel::SerialParam(),
  dir = "~/livemount/Bio_data/ORFik_experiments"
)
valid_collections <- c("all_samples-Homo_sapiens")
metadata <- fread("~/livemount/Bio_data/NGS_pipeline/metadata_rc.csv")
all_exp_meta <- all_exp[
  name %in% c(grep("all_samples-", name), valid_collections),
]
RiboCrypt_app(
  FALSE,
  all_exp = all_exp_meta,
  metadata = metadata,
  browser_options = c(
    default_experiment = "all_samples-Homo_sapiens",
    default_gene = "ATF4-ENSG00000128272",
    default_gene_meta = "ATF4-ENSG00000128272",
    default_kmer = 9,
    default_frame_type = "area",
    plot_on_start = FALSE
  )
)
