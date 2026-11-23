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

test_that("RGST status gets a colour and legend entry; uncoloured statuses still render", {
  meta <- list(status_levels = c("understood", "RGST", "not understood", "new-status"),
               status_colors = c("understood" = "#a1d99b", "not understood" = "#de2d26"))
  cols <- RiboCrypt:::explain_status_colors(meta)
  expect_true(all(meta$status_levels %in% names(cols)))
  expect_equal(unname(cols["RGST"]), "#6baed6")
  expect_equal(unname(cols["new-status"]), "#bdbdbd")
  expect_equal(names(cols)[1:4], meta$status_levels)
  lab <- RiboCrypt:::explain_status_labels(meta)
  expect_match(lab[["RGST"]], "reproducible gene-specific trend")
  meta$status_labels <- c("RGST" = "RGST (renamed)")
  expect_equal(RiboCrypt:::explain_status_labels(meta)[["RGST"]], "RGST (renamed)")
  html <- as.character(RiboCrypt:::explain_status_legend_ui(meta))
  expect_match(html, "RGST \\(renamed\\)"); expect_match(html, "#6baed6"); expect_match(html, "new-status")
  dt <- data.table::data.table(x0 = c(1, 11, 21), x1 = c(10, 20, 30))
  tr <- RiboCrypt:::explain_bar_traces(dt, c("RGST", "new-status", "understood"), cols, rep(0.9, 3), rep(0, 3), c("a", "b", "c"))
  expect_setequal(vapply(tr, `[[`, "", "name"), c("RGST", "new-status", "understood"))
})

# ---- translon calls lane ----------------------------------------------------

write_min_store <- function(dir) {
  dir.create(dir, showWarnings = FALSE)
  fst::write_fst(data.frame(transcript = "TX1", gene = "G1", n = 200L, cs = 11L, ce = 190L,
                            w_from = 1L, w_to = 1L, l_from = 1L, l_to = 0L, p_from = 1L, p_to = 0L,
                            f_from = 1L, f_to = 0L), file.path(dir, "index.fst"))
  fst::write_fst(data.frame(start = 1L, end = 200L, region = 2L, status = 1L, type = NA_integer_, O = 10, E = 10,
                            usable_frac = 1, status_all = 1L, type_all = NA_integer_, O_all = 10, E_all = 10),
                 file.path(dir, "windows.fst"))
  fst::write_fst(data.frame(start = integer(), end = integer()), file.path(dir, "lowmap.fst"))
  fst::write_fst(data.frame(codon = integer(), plddt = integer()), file.path(dir, "plddt.fst"))
  fst::write_fst(data.frame(s = integer(), e = integer(), ftype = integer(), note = character()),
                 file.path(dir, "features.fst"))
  saveRDS(list(organism = "Homo sapiens", status_levels = "understood", status_colors = c(understood = "#a1d99b"),
               region_levels = c("leader", "CDS", "trailer"), type_levels = character(),
               ftype_levels = "domain", feature_lanes = list("Domains" = "domain")), file.path(dir, "meta.rds"))
  dir
}

write_translons <- function(dir, ver = "v4_3") {
  # TX1: uORF 21-29 (high) and overlapping uoORF 26-120 crossing the exon junction at 100 (candidate, slot 2),
  # CTE 191-200 (high). TX2 (not in the main store): intORF 5-40.
  calls <- data.frame(start = c(21L, 26L, 191L, 5L), end = c(29L, 120L, 200L, 40L), slot = c(1L, 2L, 1L, 1L),
                      type = c(1L, 2L, 3L, 4L), tier = c(1L, 2L, 1L, 2L), start_codon = c(2L, 1L, 3L, 1L),
                      P = c(0.99, 0.6, NA, 0.7), q = c(0.01, 0.15, 0, 0.12), conserved = c(TRUE, FALSE, NA, FALSE),
                      rt_eff = c(NA, NA, 0.0465, NA), id = c("a", "b", "c", "d"),
                      genomic = c("1:1021-1029:-", "1:1401-1420,1001-1080:-", "1:1401-1410:-", "1:2005-2040:+"))
  fst::write_fst(calls, file.path(dir, paste0("translons_", ver, ".fst")))
  fst::write_fst(data.frame(transcript = c("TX1", "TX2"), gene = c("G1", "G2"), n = c(200L, 100L), chr = "1",
                            strand = c("-", "+"), from = c(1L, 4L), to = c(3L, 4L), n_slots = c(2L, 1L),
                            n_high = c(2L, 0L)), file.path(dir, paste0("translons_", ver, "_index.fst")))
  saveRDS(list(version = ver, type_levels = c("uORF", "uoORF", "CTE", "intORF"),
               type_colors = c(uORF = "#E69F00", uoORF = "#D55E00", CTE = "#CC79A7", intORF = "#0072B2"),
               tier_levels = c("high", "candidate"), tier_opacity = c(high = 0.95, candidate = 0.4),
               codon_levels = c("ATG", "CTG", "none")), file.path(dir, paste0("translons_", ver, "_meta.rds")))
}

test_that("stores without translon files work as before (no lane, no layer choice)", {
  dir <- write_min_store(tempfile("explain_store_old"))
  store <- RiboCrypt:::explain_store_open(dir)
  expect_null(store$translons)
  expect_false("translons" %in% RiboCrypt:::explain_layer_choices(store))
  expect_false("translons" %in% RiboCrypt:::explain_layer_choices())
  tx <- GRanges("chr1", IRanges(c(1401, 1001), c(1500, 1100)), "-")
  explain <- list(active = TRUE, ids = "TX1", layers = c("status", "translons"), tx_exons = GRangesList(TX1 = tx))
  ex <- RiboCrypt:::explain_locus_data(explain, GRangesList(TX1 = genomic_display(tx)), store)
  expect_equal(nrow(ex$translons), 0)
  expect_named(RiboCrypt:::explain_bottom_panels(ex, GRangesList(TX1 = genomic_display(tx))), "explain_status")
})

test_that("translon lane: discovery, lazy per-transcript reads, junction split, hover and tier opacity", {
  dir <- write_min_store(tempfile("explain_store_tl"))
  write_translons(dir, "v4_2")
  write_translons(dir, "v4_10")  # newest version wins (numeric, not lexical)
  store <- RiboCrypt:::explain_store_open(dir)
  expect_equal(store$translons$version, "v4_10")
  expect_true("translons" %in% RiboCrypt:::explain_layer_choices(store))
  expect_equal(nrow(RiboCrypt:::explain_translons_tx(store$translons, "TX1")$calls), 3)
  expect_null(RiboCrypt:::explain_translons_tx(store$translons, "nope"))
  # translon-only transcripts become candidates only with the layer selected
  expect_setequal(RiboCrypt:::explain_known_ids(store, "status"), "TX1")
  expect_setequal(RiboCrypt:::explain_known_ids(store, c("status", "translons")), c("TX1", "TX2"))
  tx1 <- GRanges("chr1", IRanges(c(1401, 1001), c(1500, 1100)), "-")
  tx2 <- GRanges("chr1", IRanges(2001, 2100), "+")
  grl <- GRangesList(TX1 = tx1, TX2 = tx2)
  expect_equal(RiboCrypt:::explain_candidate_ids(store, c("TX2"), grl, GRangesList(), layers = "translons"), "TX2")
  expect_length(RiboCrypt:::explain_candidate_ids(store, c("TX2"), grl, GRangesList(), layers = "status"), 0)
  explain <- list(active = TRUE, ids = c("TX1"), layers = c("status", "translons"), tx_exons = grl["TX1"])
  display <- GRangesList(TX1 = genomic_display(tx1))
  ex <- RiboCrypt:::explain_locus_data(explain, display, store)
  tl <- ex$translons
  # uoORF 26-120 crosses the junction after tx nt 100 -> two boxes, same id
  expect_equal(nrow(tl[id == "b"]), 2)
  expect_equal(tl[id == "b"]$x0, c(26, 401))
  expect_equal(tl[id == "b"]$x1, c(100, 420))
  ref <- reference_x(26, 120, tx1, genomic_display(tx1))[[1]]
  expect_equal(unname(as.matrix(tl[id == "b", .(x0, x1)])), unname(ref))
  panels <- RiboCrypt:::explain_bottom_panels(ex, display)
  expect_true("explain_translons" %in% names(panels))
  tr <- panels$explain_translons$x$data
  expect_setequal(vapply(tr, `[[`, "", "name"), c("uORF", "uoORF", "CTE"))
  cte <- tr[[which(vapply(tr, `[[`, "", "name") == "CTE")]]
  expect_match(cte$hovertext, "readthrough efficiency: 0.0465")
  expect_match(cte$hovertext, "CTE</b> \\(high tier\\)")
  expect_match(cte$hovertext, "start codon none")
  uo <- tr[[which(vapply(tr, `[[`, "", "name") == "uoORF")]]
  expect_equal(as.numeric(uo$marker$opacity), c(0.4, 0.4))
  expect_match(uo$hovertext[1], "conserved translon \\(mouse\\): no")
  expect_match(uo$hovertext[1], "1:1401-1420,1001-1080:-")
  u <- tr[[which(vapply(tr, `[[`, "", "name") == "uORF")]]
  expect_equal(as.numeric(u$marker$opacity), 0.95)
  expect_match(u$hovertext, "q = 0.01, P = 0.99")
  # slot 1 drawn above slot 2 within the transcript block (2 rows)
  expect_gt(as.numeric(u$base), as.numeric(uo$base[1]))
  expect_equal(panels$explain_translons$x$layout$yaxis$range, c(0, 2))
  # a translon-only transcript (TX2, not in the main store) still gets its lane; store lanes unaffected
  explain2 <- list(active = TRUE, ids = c("TX2", "TX1"), layers = c("status", "translons"), tx_exons = grl)
  display2 <- GRangesList(TX2 = GRanges("chr1", IRanges(1001, 2100), "+"))
  ex2 <- RiboCrypt:::explain_locus_data(explain2, display2, store)
  expect_setequal(unique(ex2$translons$tx), c("TX1", "TX2"))
  expect_equal(unique(ex2$status$lane), 1L)
  expect_equal(ex2$translons[tx == "TX2"]$x0, 1005)
  p2 <- RiboCrypt:::explain_bottom_panels(ex2, display2)$explain_translons
  expect_equal(p2$x$layout$yaxis$range, c(0, 3))
  # length mismatch with the annotation: lane skipped with a warning
  bad <- GRangesList(TX2 = GRanges("chr1", IRanges(2001, 2050), "+"))
  expect_warning(res <- RiboCrypt:::explain_translons_locus(store$translons, "TX2", bad[[1]], bad[[1]]), "length mismatch")
  expect_null(res)
  html <- as.character(RiboCrypt:::explain_translons_legend_ui(store$translons$meta))
  expect_match(html, "v4_10"); expect_match(html, "#CC79A7")
})
