# Source with devtools::load_all("."); real ATF4, T + TC predictions, full-library clusters.
local({
  df <- ORFik::read.experiment("all_samples-Homo_sapiens", validate = FALSE)
  metadata <- data.table::fread("/media/roler/S/data/Bio_data/projects/metadata_done_samples_extended_qc.csv")
  genes <- RiboCrypt:::get_gene_name_categories(df)
  id <- RiboCrypt:::tx_from_gene_list(genes, "ATF4-ENSG00000128272")[1]
  annotation <- RiboCrypt:::subset_tx_by_region(df, id, "leader+cds", ORFik::loadRegion(df, "cds"), ORFik::loadRegion(df))
  controller <- list(dff = df, id = id, display_region = annotation$region, tx_annotation = annotation$region,
    annotation = RiboCrypt:::observed_cds_annotation_internal(id, annotation$cds_annotation, FALSE),
    customRegions = RiboCrypt:::observed_translon_annotation(id, df, FALSE, TRUE, TRUE),
    viewMode = FALSE, collapsed_introns_width = 0,
    table_path = RiboCrypt:::collection_path_from_exp(df, id, grl_all = annotation$region))
  set.seed(42)
  full <- RiboCrypt:::compute_collection_table(controller$table_path, RiboCrypt:::get_lib_sizes_file(df), df,
    c("TISSUE", "CELL_LINE"), "maxNormalized", 1, metadata, min_count = 100,
    as_list = TRUE, enrichment_term = "Clusters", clusters = 5)
  grouped <- RiboCrypt:::allsamples_metadata_clustering(full, "Clusters", compute_stats = FALSE)
  time <- system.time(result <- RiboCrypt:::megabrowser_translon_analysis(controller, full, grouped))
  stopifnot(any(result$regions$Type == "clean_cds"), any(grepl("T[0-9].*TC[0-9]", result$regions$Labels)),
    all(result$statistics$Libraries[1:5] > 0), all(result$ratios$Run %in% colnames(full$table)))
  saveRDS(list(controller = controller, result = result, full = full, grouped = grouped, elapsed = time),
    "/tmp/megabrowser-atf4-translons.rds")
  print(id)
  print(result$regions)
  print(time)
  print(head(result$statistics))
})
