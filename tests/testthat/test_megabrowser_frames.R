test_that("frame bins retain separate phase sums and conserve counts", {
  mat <- matrix(seq_len(40), nrow = 20, dimnames = list(NULL, c("a", "b")))
  for (k in c(2L, 3L, 4L, 9L, 10L)) {
    result <- RiboCrypt:::multiSampleBinRows(mat, k, TRUE)
    groups <- ((seq_len(nrow(mat)) - 1L) %/% (3L * ceiling(k / 3))) * 3L +
      ((seq_len(nrow(mat)) - 1L) %% 3L) + 1L
    expect_equal(unname(result), unname(rowsum(mat, groups, reorder = TRUE)))
    expect_equal(unname(colSums(result)), unname(colSums(mat)))
  }
  expect_identical(RiboCrypt:::multiSampleBinRows(mat, 30L, TRUE), mat)
  expect_equal(unname(RiboCrypt:::multiSampleBinRows(mat, 9L, FALSE)),
               unname(rowsum(mat, (seq_len(20) - 1L) %/% 9L + 1L)))
  expect_equal(nrow(RiboCrypt:::multiSampleBinRows(matrix(1, 180, 1), 9L, TRUE)), 60L)
  expect_equal(nrow(RiboCrypt:::multiSampleBinRows(matrix(1, 180, 1), 9L, FALSE)), 20L)
})

test_that("translon and cluster collapse retain three distinct phase cells", {
  raw <- matrix(seq_len(48), nrow = 12, dimnames = list(NULL, letters[1:4]))
  table <- matrix(1, 4, 4, dimnames = list(NULL, letters[1:4]))
  attr(table, "split_by_frame") <- TRUE
  attr(table, "ratio") <- 3L
  attr(table, "km") <- list(cluster = c(1L, 1L, 2L, 2L))
  attr(table, "row_order_list") <- list(c(1L, 2L), c(3L, 4L))
  attr(table, "clusters") <- 2L
  regions <- list(R1 = list(labels = "uORF", kind = "Translon", ranges = IRanges::IRanges(1, 3)),
                  R2 = list(labels = "CDS", kind = "CDS", ranges = IRanges::IRanges(5, 12)))
  view <- RiboCrypt:::megabrowser_translon_display(table, raw, regions)
  expect_identical(rownames(view), c("R1_F0", "R1_F1", "R1_F2", "R2_F0", "R2_F1", "R2_F2"))
  expect_equal(as.vector(view["R2_F0", ]), as.vector(colMeans(raw[c(5, 8, 11), ])))
  expect_equal(attr(view, "translon_regions")$Length_nt, c(1, 1, 1, 3, 3, 2))
  meta <- data.table::data.table(Run = letters[1:4], cluster = c(1, 1, 2, 2), grouping = "x")
  display <- RiboCrypt:::megabrowser_collapsed_display(view, meta)
  groups <- RiboCrypt:::megabrowser_display_groups(meta, colnames(view))
  for (g in names(groups)) expect_equal(as.vector(display$table[, g]),
    as.vector(rowMeans(view[, groups[[g]], drop = FALSE])))
  scored <- RiboCrypt:::megabrowser_translon_cluster_score(display)
  d <- RiboCrypt:::megabrowser_plotly_heatmap_data(scored$table)
  expect_true(d$frame_view)
  expect_equal(unname(floor(d$z[1, ])), rep(0:2, 2))
  expect_match(d$hovertemplate, "Mean coverage")
  expect_match(d$text[1, 1], "frame 0")
  expect_equal(nrow(d$z), 2)
  expect_equal(ncol(d$z), 6)
  expect_match(RiboCrypt:::megabrowser_display_counts(table, scored), "2 regions \\(6 phase cells\\)")
  click <- RiboCrypt:::megabrowser_cell_selection(list(x = 5, y = 1), view, scored, meta)
  expect_identical(click$region, "R2_F1")
  details <- RiboCrypt:::megabrowser_cell_libraries(click, view, meta, data.table::data.table(Run = letters[1:4]))
  expect_equal(details$Coverage_sum, as.vector(view["R2_F1", ]) * 3)
  expect_s4_class(RiboCrypt:::megabrowser_complex_heatmap(scored$table, c("white", "blue")), "Heatmap")
})

test_that("clean CDS keeps its parent frame through discontiguous intervals", {
  raw <- matrix(seq_len(24), nrow = 12, dimnames = list(NULL, c("a", "b")))
  table <- raw
  attr(table, "split_by_frame") <- TRUE
  attr(table, "km") <- list(cluster = c(1, 1))
  attr(table, "row_order_list") <- list(1:2)
  regions <- list(R1 = list(labels = "CDS", kind = "CDS", ranges = IRanges::IRanges(2, 12)),
    R2 = list(labels = "clean CDS", kind = "clean_cds", parent_labels = "CDS",
      parent_start = 2, parent_end = 12, ranges = IRanges::IRanges(c(4, 9), c(5, 12))))
  view <- RiboCrypt:::megabrowser_translon_display(table, raw, regions)
  expect_identical(rownames(view), c("R2_F0", "R2_F1", "R2_F2"))
  expect_equal(attr(view, "translon_regions")$Frame_origin, rep(2, 3))
  expect_equal(as.vector(view["R2_F0", ]), as.vector(colMeans(raw[c(5, 11), ])))
  expect_equal(as.vector(view["R2_F1", ]), as.vector(colMeans(raw[c(9, 12), ])))
  expect_equal(as.vector(view["R2_F2", ]), as.vector(colMeans(raw[c(4, 10), ])))
})

test_that("short regions leave absent phases undefined rather than inventing zero coverage", {
  raw <- matrix(seq_len(12), nrow = 6, dimnames = list(NULL, c("a", "b")))
  attr(raw, "split_by_frame") <- TRUE
  attr(raw, "km") <- list(cluster = c(1, 1))
  attr(raw, "row_order_list") <- list(1:2)
  regions <- list(R1 = list(labels = "short", kind = "Translon", ranges = IRanges::IRanges(1, 2)))
  view <- RiboCrypt:::megabrowser_translon_display(raw, raw, regions)
  expect_equal(attr(view, "translon_regions")$Length_nt, c(1, 1, 0))
  expect_true(all(is.na(view["R1_F2", ])))
  d <- RiboCrypt:::megabrowser_plotly_heatmap_data(view)
  expect_true(all(is.na(d$z[, 3])))
})

test_that("frame normalization retains its phase and coordinate metadata", {
  mat <- matrix(seq_len(40), nrow = 20, dimnames = list(NULL, c("a", "b")))
  split <- RiboCrypt:::normalize_collection(mat, "transcriptNormalized", kmer = 9, split_by_frame = TRUE)
  expect_true(attr(split, "split_by_frame"))
  expect_equal(attr(split, "ratio"), 3)
  expect_equal(nrow(split), 8L)
  expect_equal(colSums(expm1(split)), c(a = 1e6, b = 1e6), tolerance = 1e-8)
  short <- RiboCrypt:::normalize_collection(mat[1:4, ], "transcriptNormalized", kmer = 9, split_by_frame = TRUE)
  expect_equal(attr(short, "ratio"), 1)
  ordinary <- RiboCrypt:::normalize_collection(mat, "transcriptNormalized", kmer = 9, split_by_frame = FALSE)
  expect_false(attr(ordinary, "split_by_frame"))
  expect_equal(attr(ordinary, "ratio"), 9)
})

test_that("frame heatmap encodes colour while keeping exact hover values", {
  mat <- matrix(c(0, 1, 2, 3, 4, 5, 6, 7), nrow = 4,
                dimnames = list(NULL, c("a", "b")))
  attr(mat, "km") <- list(cluster = c(1L, 1L))
  attr(mat, "row_order_list") <- list(c(2L, 1L))
  attr(mat, "ratio") <- 1L
  attr(mat, "split_by_frame") <- TRUE
  d <- RiboCrypt:::megabrowser_plotly_heatmap_data(mat)
  expect_true(d$frame_view)
  expect_equal(unname(d$customdata), unname(t(mat)[c(2L, 1L), ]))
  expect_equal(floor(d$z[1, ]), c(0, 1, 2, 0))
  expect_equal(d$x, 1:4)
  expect_match(d$hovertemplate, "customdata", fixed = TRUE)
  uniform <- RiboCrypt:::megabrowser_frame_encode(matrix(5, 2, 3))
  expect_gt(uniform[1, 1], 0.99)
  for (template in list(NULL, RiboCrypt:::megabrowserHeatmapPlotlyTemplate())) {
    p <- plotly::plotly_build(RiboCrypt:::megabrowser_plotly_heatmap(mat, c("white", "blue"), template))
    expect_equal(p$x$data[[1]]$zmin, 0)
    expect_equal(p$x$data[[1]]$zmax, 3)
    expect_equal(as.vector(p$x$data[[1]]$customdata), as.vector(d$customdata))
    expect_equal(as.vector(p$x$data[[1]]$z), as.vector(d$z))
  }
  attr(mat, "split_by_frame") <- FALSE
  attr(mat, "collapsed_translons") <- TRUE
  attr(mat, "translon_regions") <- data.frame(Labels = rep("region", nrow(mat)))
  expect_false(RiboCrypt:::megabrowser_plotly_heatmap_data(mat)$frame_view)
})

test_that("CDS phase follows spliced coding sequence on both strands", {
  gr <- function(start, end, strand) GenomicRanges::GRanges("chr1", IRanges::IRanges(start, end), strand)
  for (strand in c("+", "-")) {
    if (strand == "+") {
      display <- GenomicRanges::GRangesList(tx = gr(c(100, 200), c(108, 207), strand))
      cds <- GenomicRanges::GRangesList(tx = gr(c(103, 200), c(108, 205), strand))
      genomic <- GenomicRanges::GRangesList(tx = gr(100, 207, strand))
      positions <- c(4:9, 10:15)
      genomic_positions <- c(4:9, 101:106)
    } else {
      display <- GenomicRanges::GRangesList(tx = gr(c(200, 100), c(208, 107), strand))
      cds <- GenomicRanges::GRangesList(tx = gr(c(200, 102), c(205, 107), strand))
      genomic <- GenomicRanges::GRangesList(tx = gr(100, 208, strand))
      positions <- c(4:9, 10:15)
      genomic_positions <- c(4:9, 102:107)
    }
    phase <- RiboCrypt:::megabrowser_cds_phase(display, cds)
    expect_equal(phase[positions], rep(0:2, 4))
    phase_genomic <- RiboCrypt:::megabrowser_cds_phase(genomic, cds)
    expect_equal(phase_genomic[genomic_positions], rep(0:2, 4))
    expect_null(RiboCrypt:::megabrowser_cds_phase(display, GenomicRanges::GRangesList(tx = gr(500, 510, strand))))
  }
})

test_that("CDS anchored bins preserve phase and counts, including empty slots", {
  phase <- c(2L, 0L, 1L, 2L, 0L, 1L, 2L, 0L, 1L, 0L)
  mat <- matrix(seq_along(phase), ncol = 1, dimnames = list(NULL, "a"))
  binned <- RiboCrypt:::multiSampleBinRows(mat, 3L, TRUE, phase)
  expect_equal(as.vector(colSums(binned)), sum(mat))
  expect_false(tail(attr(binned, "frame_cell_present"), 1))
  normalized <- RiboCrypt:::normalize_collection(mat, "transcriptNormalized", kmer = 3, split_by_frame = TRUE, frame_phase = phase)
  expect_identical(attr(normalized, "frame_reference"), "CDS")
  expect_equal(attr(normalized, "frame_phase"), rep(0:2, 4))
  expect_equal(as.vector(sum(expm1(normalized))), 1e6, tolerance = 1e-8)
  raw <- RiboCrypt:::normalize_collection(mat, "transcriptNormalized", kmer = 1, split_by_frame = TRUE, frame_phase = phase)
  encoded <- RiboCrypt:::megabrowser_frame_encode(t(raw), frame_phase = attr(raw, "frame_phase"))
  expect_equal(as.integer(floor(encoded)), phase)
})

test_that("a clipped view uses the original CDS frame when the TIS is outside it", {
  gr <- function(s, e, strand) GenomicRanges::GRanges("chr1", IRanges::IRanges(s, e), strand)
  for (strand in c("+", "-")) {
    cds <- GenomicRanges::GRangesList(tx = gr(100, 111, strand))
    display <- GenomicRanges::GRangesList(tx = if (strand == "+") gr(104, 109, strand) else gr(102, 107, strand))
    expect_identical(RiboCrypt:::megabrowser_cds_phase(display, cds), c(1L, 2L, 0L, 1L, 2L, 0L))
  }
})
