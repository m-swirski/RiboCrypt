library(GenomicRanges)
library(data.table)

# Two-exon transcripts: exon 1 = 100 nt, exon 2 = 100 nt, intron 300 nt
plus_tx <- GRanges("chr1", IRanges(c(1001, 1401), c(1100, 1500)), "+")
minus_tx <- GRanges("chr1", IRanges(c(1401, 1001), c(1500, 1100)), "-")  # ORFik order: 5' exon first
genomic_display <- function(tx) GRanges("chr1", IRanges(min(start(tx)), max(end(tx))), strand(tx)[1])

# Reference: ORFik maps transcript ranges to genome, then genome to display x
reference_x <- function(start, end, tx, display) {
  grl <- GRangesList(tx)
  g <- ORFik::pmapFromTranscriptF(IRanges(start, end), grl[rep(1, length(start))], removeEmpty = TRUE)
  d <- display
  lapply(seq_along(g), function(i) {
    gi <- g[[i]]
    x <- if (as.character(strand(d))[1] == "-") {
      cbind(end(d) - end(gi) + 1, end(d) - start(gi) + 1)
    } else cbind(start(gi) - start(d) + 1, end(gi) - start(d) + 1)
    x[order(x[, 1]), , drop = FALSE]
  })
}

test_that("transcript view maps windows to themselves (plus and minus strand)", {
  for (tx in list(plus_tx, minus_tx)) {
    map <- RiboCrypt:::explain_tx_display_map(tx, tx)
    res <- RiboCrypt:::explain_map_ranges(c(1, 91, 181), c(90, 180, 200), map)
    expect_equal(res$id, 1:3)
    expect_equal(res$x0, c(1, 91, 181))
    expect_equal(res$x1, c(90, 180, 200))
  }
})

test_that("genomic view splits a junction-crossing window (plus strand)", {
  d <- genomic_display(plus_tx)
  map <- RiboCrypt:::explain_tx_display_map(plus_tx, d)
  res <- RiboCrypt:::explain_map_ranges(91, 180, map)  # crosses the junction after tx nt 100
  expect_equal(nrow(res), 2)
  expect_equal(res$x0, c(91, 401))
  expect_equal(res$x1, c(100, 480))
  ref <- reference_x(91, 180, plus_tx, d)[[1]]
  expect_equal(unname(as.matrix(res[, .(x0, x1)])), unname(ref))
})

test_that("genomic view splits a junction-crossing window (minus strand)", {
  d <- genomic_display(minus_tx)
  map <- RiboCrypt:::explain_tx_display_map(minus_tx, d)
  res <- RiboCrypt:::explain_map_ranges(c(1, 91), c(90, 180), map)
  ref <- reference_x(c(1, 91), c(90, 180), minus_tx, d)
  expect_equal(unname(as.matrix(res[id == 1, .(x0, x1)])), unname(ref[[1]]))
  expect_equal(unname(as.matrix(res[id == 2, .(x0, x1)])), unname(ref[[2]]))
  # display runs 5'->3' of the minus-strand gene: tx nt 1 is display x 1
  expect_equal(res[id == 1]$x0, 1)
  expect_equal(res[id == 2]$x0, c(91, 401))
  expect_equal(res[id == 2]$x1, c(100, 480))
})

test_that("collapsed introns and partial displays map consistently", {
  # pseudo-intron display: exons with 20 nt flanks into the intron
  d <- GRanges("chr1", IRanges(c(1001, 1381), c(1120, 1500)), "+")
  map <- RiboCrypt:::explain_tx_display_map(plus_tx, d)
  res <- RiboCrypt:::explain_map_ranges(91, 180, map)
  expect_equal(res$x0, c(91, 141))
  expect_equal(res$x1, c(100, 220))
  # display covering only exon 2: window part on exon 1 dropped
  d2 <- GRanges("chr1", IRanges(1401, 1500), "+")
  map2 <- RiboCrypt:::explain_tx_display_map(plus_tx, d2)
  res2 <- RiboCrypt:::explain_map_ranges(c(1, 91), c(90, 180), map2)
  expect_equal(res2$id, 2L)
  expect_equal(c(res2$x0, res2$x1), c(1, 80))
})

test_that("point mapping agrees with range mapping on both strands", {
  for (tx in list(plus_tx, minus_tx)) {
    d <- genomic_display(tx)
    map <- RiboCrypt:::explain_tx_display_map(tx, d)
    pos <- c(1, 50, 100, 101, 200)
    x <- RiboCrypt:::explain_map_points(pos, map)
    r <- RiboCrypt:::explain_map_ranges(pos, pos, map)
    expect_equal(x, r$x0)
    expect_true(is.na(RiboCrypt:::explain_map_points(250, map)))
  }
})

test_that("expected line is scaled to the library coverage at usable positions", {
  expected <- data.table(x0 = c(1L, 11L), x1 = c(10L, 20L), density = c(1, 3),
                         window = 1:2, E_window = c(10, 30), E_shown_total = 40)
  profile <- data.table(position = 1:20, count = rep(c(2, 6), each = 10), frame = factor(0))
  line <- RiboCrypt:::explain_expected_line(expected, profile)
  # observed 80 vs E 40 -> scale 2
  expect_equal(line$y[c(1, 4)], c(2, 6))
  expect_equal(line$x[1:2], c(0.5, 10.5))
  expect_true(is.na(line$x[3]))
  log_line <- RiboCrypt:::explain_expected_line(expected, copy(profile)[, count := log2(count + 1)], TRUE)
  expect_equal(log_line$y[1], log2(3))
})

test_that("store layers render as panels for a synthetic transcript", {
  dir <- tempfile("explain_store")
  dir.create(dir)
  fst::write_fst(data.frame(transcript = "TX1", gene = "G1", n = 200L, cs = 11L, ce = 190L,
                            w_from = 1L, w_to = 3L, l_from = 1L, l_to = 1L, p_from = 1L, p_to = 3L,
                            f_from = 1L, f_to = 1L), file.path(dir, "index.fst"))
  fst::write_fst(data.frame(start = c(1L, 11L, 101L), end = c(10L, 100L, 200L), region = 1:3,
                            status = c(1L, 2L, 5L), type = c(NA, 1L, NA), O = c(5, 500, 50), E = c(5, 250, 60),
                            usable_frac = c(1, 0.9, 1), status_all = c(1L, 2L, 1L), type_all = NA_integer_,
                            O_all = c(5, 520, 55), E_all = c(5, 260, 60)), file.path(dir, "windows.fst"))
  fst::write_fst(data.frame(start = 50L, end = 58L), file.path(dir, "lowmap.fst"))
  fst::write_fst(data.frame(codon = 1:3, plddt = c(950L, 600L, 300L)), file.path(dir, "plddt.fst"))
  fst::write_fst(data.frame(s = 1L, e = 2L, ftype = 1L, note = "test"), file.path(dir, "features.fst"))
  saveRDS(list(organism = "Homo sapiens", status_levels = c("understood", "not understood", "study-specific",
                 "family-level", "unknown", "contamination"),
               status_colors = c("understood" = "#a1d99b", "not understood" = "#de2d26",
                                 "study-specific" = "#fd8d3c", "family-level" = "#9e9ac8",
                                 "unknown" = "#d9d9d9", "contamination" = "#8c510a"),
               region_levels = c("leader", "CDS", "trailer"), type_levels = "CDS density bump",
               ftype_levels = "domain", feature_lanes = list("Domains" = "domain")),
          file.path(dir, "meta.rds"))
  store <- RiboCrypt:::explain_store_open(dir)
  expect_true(RiboCrypt:::explain_organism_supported("Homo_sapiens", store))
  expect_false(RiboCrypt:::explain_organism_supported("Mus musculus", store))
  tx <- GRanges("chr1", IRanges(c(1001, 1401), c(1100, 1500)), "-")
  tx <- tx[order(start(tx), decreasing = TRUE)]
  explain <- list(active = TRUE, ids = "TX1", layers = c("status", "expected", "plddt", "features"),
                  tx_exons = GRangesList(TX1 = tx))
  display <- GRangesList(TX1 = genomic_display(tx))
  ex <- RiboCrypt:::explain_locus_data(explain, display, store)
  # the CDS window 11-100 does not cross the junction (exon 1 = tx 1-100); the trailer 101-200 is exon 2
  expect_equal(nrow(ex$status), 3)
  expect_equal(ex$status[start == 101]$x0, 401)
  # expected segments skip the lowmap interval 50-58
  expect_true(all(ex$expected$x1 < 50 | ex$expected$x0 > 58))
  expect_equal(ex$expected[window == 2]$density[1], 250 / (90 * 0.9))
  panels <- RiboCrypt:::explain_bottom_panels(ex, display)
  expect_named(panels, c("explain_status", "explain_plddt", "explain_features"))
  built <- panels$explain_status
  expect_true(any(grepl("log2 O/E", unlist(lapply(built$x$data, `[[`, "hovertext")))))
  expect_equal(RiboCrypt:::explain_store_tx(store, "TX1")$tx$gene, "G1")
  expect_null(RiboCrypt:::explain_store_tx(store, "missing"))
})
