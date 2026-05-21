make_exon_coverage_fixture <- function() {
  cds <- GenomicRanges::GRangesList(
    txA = GenomicRanges::GRanges(
      "chr1",
      IRanges::IRanges(c(100, 200), c(102, 203)),
      "+"
    ),
    txB = GenomicRanges::GRanges(
      "chr1",
      IRanges::IRanges(300, 303),
      "+"
    )
  )
  high_exon_reads <- rep(100:102, 5)
  reads <- GenomicRanges::GRanges(
    "chr1",
    IRanges::IRanges(c(high_exon_reads, 201, 301), width = 1),
    "+"
  )
  list(cds = cds, reads = reads)
}

test_that("coveragePerExonTiling keeps exon ids after transcript tiling", {
  fixture <- make_exon_coverage_fixture()

  tiled <- coveragePerExonTiling(fixture$cds, fixture$reads)

  expect_true(all(c(
    "transcript_id", "exon_id", "exon_rank", "exon_position"
  ) %in% names(tiled)))
  expect_identical(unique(tiled[transcript_id == "txA", exon_id]),
                   c("txA:exon1", "txA:exon2"))
  expect_identical(tiled[transcript_id == "txA" & position == 4,
                         exon_position][[1]], 1L)
  expect_equal(tiled[transcript_id == "txB" & position == 2,
                     count][[1]], 1)
})

test_that("cdsExonCoverageQC flags low aggregate CDS exon coverage", {
  fixture <- make_exon_coverage_fixture()

  qc <- cdsExonCoverageQC(
    fixture$cds,
    fixture$reads,
    low_exon_fraction = 0.1,
    min_transcript_counts = 10
  )

  expect_true(qc[transcript_id == "txA" & exon_rank == 2,
                 low_coverage_exon][[1]])
  expect_true(qc[transcript_id == "txA",
                 unique(transcript_has_low_coverage_exon)])
  expect_false(qc[transcript_id == "txB",
                  exon_qc_measurable][[1]])
  expect_false(qc[transcript_id == "txB",
                  transcript_has_low_coverage_exon][[1]])
})
