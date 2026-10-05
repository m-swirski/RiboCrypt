# Run from the repo with devtools::load_all() before sourcing this file.
local({
  df <- ORFik::read.experiment("all_samples-Homo_sapiens", validate = FALSE)
  metadata <- data.table::fread("/media/roler/S/data/Bio_data/projects/metadata_done_samples_extended_qc.csv")
  genes <- RiboCrypt:::get_gene_name_categories(df)
  tx <- RiboCrypt:::tx_from_gene_list(genes, "AMD1-ENSG00000123505")[1]
  region <- RiboCrypt:::subset_tx_by_region(df, tx, "leader+cds", ORFik::loadRegion(df, "cds"), ORFik::loadRegion(df))$region
  path <- RiboCrypt:::collection_path_from_exp(df, tx, grl_all = region)
  sizes <- RiboCrypt:::get_lib_sizes_file(df)
  trials <- list()
  Rprof("/tmp/megabrowser-backend-rprof.out", interval = 0.01)
  for (trial in 1:3) {
    set.seed(42)
    load_time <- system.time(full <- RiboCrypt:::compute_collection_table(
      path, sizes, df, c("TISSUE", "CELL_LINE"), "maxNormalized", 1,
      metadata, min_count = 100, as_list = TRUE, enrichment_term = "Clusters", clusters = 5))["elapsed"]
    group_time <- system.time(grouped <- RiboCrypt:::allsamples_metadata_clustering(full, "Clusters"))["elapsed"]
    collapse_time <- system.time(display <- RiboCrypt:::megabrowser_collapsed_display(full$table, grouped$meta))["elapsed"]
    stats_time <- system.time(stats <- RiboCrypt:::megabrowser_metadata_enrichment(grouped, metadata, "CELL_LINE"))["elapsed"]
    for (collapsed in c(FALSE, TRUE)) {
      matrix <- if (collapsed) display$table else full$table
      build_time <- system.time(plot <- plotly::plotly_build(RiboCrypt:::get_meta_browser_plot(
        matrix, "default (White-Blue)", 3, "plotly")))["elapsed"]
      json_time <- system.time(json <- plotly::plotly_json(plot, jsonedit = FALSE, pretty = FALSE))["elapsed"]
      trials[[length(trials) + 1L]] <- list(trial = trial, collapsed = collapsed,
        gene = "AMD1", tx = tx, positions = nrow(matrix), libraries = ncol(full$table), display_rows = ncol(matrix),
        load_cluster_seconds = unname(load_time), grouping_seconds = unname(group_time),
        collapse_seconds = unname(collapse_time), statistics_seconds = unname(stats_time),
        plot_build_seconds = unname(build_time), json_seconds = unname(json_time), json_bytes = nchar(json, type = "bytes"))
    }
  }
  Rprof(NULL)
  jsonlite::write_json(list(R = R.version.string, trials = trials), "/tmp/megabrowser-backend.json", auto_unbox = TRUE, pretty = TRUE)
  print(trials)
})
