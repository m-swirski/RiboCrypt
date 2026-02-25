test_that("get_ratio_interval parses and validates inputs", {
  expect_equal(RiboCrypt:::get_ratio_interval("5:10;20:30"), c(5, 10, 20, 30))
  expect_equal(RiboCrypt:::get_ratio_interval("7"), c(7, 7))
  expect_error(RiboCrypt:::get_ratio_interval("5:10:15"))
  expect_error(suppressWarnings(RiboCrypt:::get_ratio_interval("a:b")))
})

test_that("multiSampleBinRows bins rows as expected", {
  mat <- matrix(1:20, nrow = 10, ncol = 2)
  binned <- RiboCrypt:::multiSampleBinRows(mat, ratio = 3)
  expect_equal(nrow(binned), 4)
  expect_equal(binned[1, ], colSums(mat[1:3, ]))

  mat_small <- matrix(1:6, nrow = 3, ncol = 2)
  expect_equal(RiboCrypt:::multiSampleBinRows(mat_small, ratio = 3), mat_small)
})

test_that("summary_track returns expected columns and values", {
  mat <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 3, ncol = 2)
  res <- RiboCrypt:::summary_track(mat)
  expect_true(all(c("count", "position", "frame") %in% names(res)))
  expect_equal(res$count, c(1 + 4, 2 + 5, 3 + 6))
  expect_equal(res$position, 1:3)
})

test_that("margin_megabrowser returns expected margins", {
  m <- RiboCrypt:::margin_megabrowser()
  expect_equal(m$l, 30)
  expect_equal(m$r, 100)
  expect_equal(m$t, 0)
  expect_equal(m$b, 0)
})

test_that("allsamples_meta_stats computes counts without clusters", {
  meta <- data.table::data.table(grouping = c("a", "a", "b"), order = 1:3, index = seq(3))
  attr(meta, "xlab") <- "X"
  attr(meta, "ylab") <- "Y"
  res <- RiboCrypt:::allsamples_meta_stats(meta)
  expect_equal(res["a", "counts"], 2)
  expect_equal(res["b", "counts"], 1)
  expect_equal(attr(res, "xlab"), "X")
  expect_equal(attr(res, "ylab"), "Y")
})

test_that("add_alpha converts colors to rgba with alpha", {
  res <- RiboCrypt:::add_alpha(c("red", "blue"), alpha = 0.5)
  expect_true(all(grepl("^rgba\\(", res)))
  expect_true(all(grepl(",0.500\\)$", res)))
})

test_that("profile_plotly_gl returns a plotly object", {
  dt <- data.table::data.table(
    position = 1:6,
    count = c(1, 2, 3, 4, 5, 6),
    frame = factor(rep(0:2, length.out = 6))
  )
  p <- RiboCrypt:::profile_plotly_gl(dt, bar_px = 2, alpha = 0.6)
  expect_true(inherits(p, "plotly"))
})

test_that("get_meta_browser_plot returns plotly heatmap for plotly type", {
  mat <- matrix(1:20, nrow = 5, ncol = 4)
  table <- data.table::data.table(mat)
  setnames(table, new = paste0("lib", seq(4)))
  data.table::setattr(table, "ratio", 1L)
  km <- stats::kmeans(t(mat), centers = 2)
  data.table::setattr(table, "km", km)
  data.table::setattr(table, "row_order_list", list("1" = seq(4)))


  p <- RiboCrypt:::get_meta_browser_plot(table, color_theme = "default (White-Blue)", plotType = "plotly")
  expect_true(inherits(p, "plotly"))

})

test_that("compute_collection_table_grouping groups metadata with fallback enrichment term", {
  metadata <- data.table::data.table(
    study_accession = c(rep("PRJNA100001", 4), rep("PRJNA100002", 4)),
    Run = c("SRR1001", "SRR1002", "SRR1003", "SRR1004", "SRR2001", "SRR2002", "SRR2003", "SRR2004"),
    BioProject = c(rep("PRJNA100001", 4), rep("PRJNA100002", 4)),
    TISSUE = c("brain", "brain", "brain", "brain", "heart", "heart", "heart", "heart"),
    CELL_LINE = c("CL1", "CL1", "CL2", "CL2", "CL3", "CL3", "CL4", "CL4"),
    CONDITION = c("ctrl", "drug", "ctrl", "drug", "ctrl", "drug", "ctrl", "drug"),
    GENE = c("GENE1", "GENE1", "GENE2", "GENE2", "GENE1", "GENE1", "GENE2", "GENE2")
  )

  table <- data.table::as.data.table(cbind(
    SRR1001 = c(130, 15, 12, 10),   # brain_ctrl
    SRR1002 = c(14, 132, 11, 10),   # brain_drug
    SRR1003 = c(128, 16, 13, 9),    # brain_ctrl
    SRR1004 = c(15, 129, 10, 11),   # brain_drug
    SRR2001 = c(12, 10, 131, 14),   # heart_ctrl
    SRR2002 = c(10, 12, 15, 130),   # heart_drug
    SRR2003 = c(13, 9, 129, 16),    # heart_ctrl
    SRR2004 = c(11, 13, 14, 128)    # heart_drug
  ))
  valid_libs <- rep(TRUE, ncol(table))
  names(valid_libs) <- colnames(table)
  data.table::setattr(table, "valid_libs", valid_libs)

  grouping_args <- list(
    metadata = metadata,
    df = NULL,
    metadata_field = c("TISSUE", "CELL_LINE", "CONDITION", "GENE"),
    table = table,
    ratio_interval = NULL,
    group_on_tx_tpm = NULL,
    decreasing_order = FALSE,
    enrichment_term = "Cluster"
  )

  res <- do.call(RiboCrypt:::compute_collection_table_grouping, grouping_args)

  expect_length(res, 8)
  expect_equal(as.character(res), c("brain", "brain", "brain", "brain", "heart", "heart", "heart", "heart"))
  expect_equal(names(res), c("SRR1001", "SRR1002", "SRR1003", "SRR1004", "SRR2001", "SRR2002", "SRR2003", "SRR2004"))
  expect_equal(attr(res, "xlab"), "TISSUE")
  expect_equal(attr(res, "meta_order"), 1:8)

  other_cols <- attr(res, "other_columns")
  expect_true(all(c("CELL_LINE", "CONDITION", "GENE", "Run") %in% colnames(other_cols)))
  expect_equal(nrow(other_cols), 8)

  run_ids <- attr(res, "runIDs")
  expect_true(all(c("Run", "BioProject") %in% colnames(run_ids)))
  expect_equal(nrow(run_ids), 8)

  # On Gene:
  grouping_args$metadata_field <- c("GENE", "TISSUE", "CELL_LINE", "CONDITION")
  res <- do.call(RiboCrypt:::compute_collection_table_grouping, grouping_args)
  expect_equal(attr(res, "meta_order"), c(1, 2, 5, 6, 3, 4, 7, 8))
  expect_equal(
    names(res),
    c("SRR1001", "SRR1002", "SRR2001", "SRR2002", "SRR1003", "SRR1004", "SRR2003", "SRR2004")
  )
  expect_equal(
    as.character(res),
    c("GENE1", "GENE1", "GENE1", "GENE1", "GENE2", "GENE2", "GENE2", "GENE2")
  )
  expect_equal(attr(res, "xlab"), "GENE")

  ## Test clustering works on the groups
  data.table::setcolorder(table, attr(res, "meta_order"))
  set.seed(1)
  table <- RiboCrypt:::clustering_megabrowser(table, clusters = 4)
  table_with_meta <- list(table = table, metadata_field = res)

  meta_and_clusters <- RiboCrypt:::allsamples_metadata_clustering(
    table_with_meta,
    grouping_args$enrichment_term
  )
  expect_true(is.list(meta_and_clusters))
  expect_true(all(c("meta", "enrich_dt") %in% names(meta_and_clusters)))
  expect_equal(nrow(meta_and_clusters$meta), 8)
  expect_true(all(c("Run", "grouping", "order", "TISSUE", "CELL_LINE", "CONDITION", "cluster", "index") %in%
                    colnames(meta_and_clusters$meta)))
  expect_equal(sort(unique(meta_and_clusters$meta$grouping)), c("GENE1", "GENE2"))
  expect_equal(attr(meta_and_clusters$meta, "xlab"), "GENE Clusters (K-means)")
  expect_equal(attr(meta_and_clusters$meta, "ylab"), "Enrichment")
  expect_true(is.data.frame(meta_and_clusters$enrich_dt))
  expect_equal(rownames(meta_and_clusters$enrich_dt), c("GENE1", "GENE2"))

  meta_tbl <- RiboCrypt:::allsamples_meta_table(meta_and_clusters)
  expect_equal(nrow(meta_tbl), 8)
  expect_equal(meta_tbl$cluster, rep(seq(4), each = 2))
  expect_true(all(c("Run", "BioProject", "grouping", "cluster", "TISSUE", "CELL_LINE", "CONDITION") %in%
                    colnames(meta_tbl)))
  expect_false(any(c("index", "order") %in% colnames(meta_tbl)))
  expect_true(any(meta_tbl$Run == "SRR1001" & meta_tbl$BioProject == "PRJNA100001" & meta_tbl$TISSUE == "brain" &
                    meta_tbl$CELL_LINE == "CL1" & meta_tbl$CONDITION == "ctrl"))

})
