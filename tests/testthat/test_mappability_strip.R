library(GenomicRanges)
library(data.table)

# Tiny store in the rc_mappability_v1 layout: chr "1", one 1 Mb bin.
write_map_store <- function(dir, exact_def, exact_nu, le_def, le_nu) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  idx <- list()
  put <- function(d, track, part) {
    d <- as.data.table(d)
    d[, chr := 1L]
    setcolorder(d, c("chr", "start", "end", "strand"))
    setorder(d, chr, start, strand)
    fst::write_fst(d, file.path(dir, paste0(track, "_", part, ".fst")))
    if (nrow(d)) idx[[length(idx) + 1L]] <<- data.table(track = track, part = part, chr = 1L, bin = 0L,
                                                         from = 1L, to = nrow(d))
  }
  put(exact_def, "exact", "def"); put(exact_nu, "exact", "nu")
  put(le_def, "le2mm", "def"); put(le_nu, "le2mm", "nu")
  fst::write_fst(rbindlist(idx), file.path(dir, "index.fst"))
  saveRDS(list(version = "rc_mappability_v1", organism = "Homo sapiens", seqlevels = c("1", "2"), bin = 1e6,
               k = 28L, psite_offset = 12L,
               tracks = list(exact = list(label = "unique (perfect match)", short = "unique",
                                          anchor = "psite_stranded", shift = c(plus = 0L, minus = 0L)),
                             le2mm = list(label = "<=2 mismatches", short = "\u22642 mm",
                                          anchor = "kmer_start_forward", shift = c(plus = 12L, minus = 15L)))),
          file.path(dir, "meta.rds"))
  dir
}

# exons 1001-1100 and 1401-1500 (transcript of 200 nt)
store_fixture <- function() {
  write_map_store(
    tempfile("rc_map"),
    exact_def = data.frame(start = c(1001L, 1401L), end = c(1100L, 1500L), strand = 1L),
    exact_nu = data.frame(start = 1051L, end = 1100L, strand = 1L, n = NA_integer_),   # 2nd half of exon 1, + strand
    le_def = data.frame(start = c(989L, 1389L), end = c(1088L, 1488L), strand = 0L),    # forward-start anchored
    le_nu = data.frame(start = c(1389L), end = c(1438L), strand = 0L, n = c(3L)))
}

test_that("genomic intervals map to display x on both strands, split at display segments", {
  d_plus <- GRanges("1", IRanges(c(1001, 1401), c(1100, 1500)), "+")
  res <- RiboCrypt:::mappability_genomic_to_x(c(1091, 1450), c(1410, 1460), d_plus)
  expect_equal(res[id == 1]$x0, c(91, 101)); expect_equal(res[id == 1]$x1, c(100, 110))
  expect_equal(c(res[id == 2]$x0, res[id == 2]$x1), c(150, 160))
  # minus strand: display runs 5'->3' of the gene, i.e. from the highest coordinate
  d_minus <- GRanges("1", IRanges(c(1401, 1001), c(1500, 1100)), "-")
  res <- RiboCrypt:::mappability_genomic_to_x(1491, 1500, d_minus)
  expect_equal(c(res$x0, res$x1), c(1, 10))
  res <- RiboCrypt:::mappability_genomic_to_x(1001, 1010, d_minus)
  expect_equal(c(res$x0, res$x1), c(191, 200))
  # genomic view with the intron shown: one display range
  res <- RiboCrypt:::mappability_genomic_to_x(1401, 1410, GRanges("1", IRanges(1001, 1500), "+"))
  expect_equal(c(res$x0, res$x1), c(401, 410))
})

test_that("block fractions: exact row strand-specific, <=2 mm row shifted to P-sites by strand", {
  store <- RiboCrypt:::mappability_store_open(store_fixture())
  d_plus <- GRanges("1", IRanges(c(1001, 1401), c(1100, 1500)), "+")
  ex <- RiboCrypt:::mappability_track_blocks(store, "exact", d_plus, max_blocks = 2L)
  # block 1 = exon 1 (half non-unique), block 2 = exon 2 (unique)
  expect_equal(ex$frac, c(0.5, 1))
  # minus-strand display: exact intervals are + strand only -> nothing defined
  d_minus <- GRanges("1", IRanges(c(1401, 1001), c(1500, 1100)), "-")
  expect_null(RiboCrypt:::mappability_track_blocks(store, "exact", d_minus))
  # <=2 mm, + strand: forward starts 989-1088 / 1389-1488 -> P-sites 1001-1100 / 1401-1500; non-unique 1401-1450
  le <- RiboCrypt:::mappability_track_blocks(store, "le2mm", d_plus, max_blocks = 2L)
  expect_equal(le$frac, c(1, 0.5))
  expect_equal(le$n_max[2], 3L)
  # minus strand: shift +15 -> defined 1004-1103 / 1404-1503 (display 1004-1100, 1404-1500), non-unique 1404-1453
  lm <- RiboCrypt:::mappability_track_blocks(store, "le2mm", d_minus, max_blocks = 200L)
  expect_equal(nrow(lm), 194)                                     # 97 + 97 defined nt, one block per nt
  x_nonuniq <- lm[frac == 0]
  expect_equal(range(x_nonuniq$x0), c(48, 97))                    # 1453..1404 on the reversed display
  # undefined positions are not drawn and not counted
  expect_true(all(lm$defined == 1))
})

test_that("strip data merges equal-class blocks; panel has two labelled rows with hover", {
  store <- RiboCrypt:::mappability_store_open(store_fixture())
  display <- GRangesList(GRanges("1", IRanges(c(1001, 1401), c(1100, 1500)), "+"))
  sd <- RiboCrypt:::mappability_strip_data(store, display, max_blocks = 200L)
  expect_named(sd$tracks, c("exact", "le2mm"))
  expect_equal(sd$tracks$exact$cls, c("100% unique", "< 20% unique", "100% unique"))
  expect_equal(sd$tracks$exact$x0, c(1, 51, 101)); expect_equal(sd$tracks$exact$x1, c(50, 100, 200))
  p <- RiboCrypt:::mappability_strip_panel(sd, display)
  expect_equal(vapply(p$x$layout$annotations, `[[`, "", "text"), c("unique", "\u22642 mm"))
  expect_equal(p$x$layout$yaxis$range, c(0, 2))
  tr <- p$x$data
  rows <- vapply(tr, function(t) t$meta$rc_row, "")
  expect_setequal(rows, c("exact", "le2mm"))
  # exact row is drawn above the <=2 mm row
  expect_gt(min(unlist(lapply(tr[rows == "exact"], `[[`, "base"))), max(unlist(lapply(tr[rows == "le2mm"], `[[`, "base"))))
  hov <- unlist(lapply(tr, `[[`, "hovertext"))
  expect_true(any(grepl("50% of P-site positions unique", hov) | grepl("0% of P-site positions unique", hov)))
  expect_true(any(grepl("max sites within 2 mm: 3", hov)))
  cols <- vapply(tr, function(t) t$marker$color, "")
  expect_true(all(cols %in% RiboCrypt:::mappability_classes()))
})

test_that("classes: light = poorly mappable, 100% only for fully unique blocks", {
  expect_equal(RiboCrypt:::mappability_class_of(c(0, 0.19, 0.5, 0.99, 1)),
               c("< 20% unique", "< 20% unique", "40-60%", "80-99%", "100% unique"))
})

test_that("chromosome naming styles are matched", {
  lv <- c("1", "X", "MT")
  expect_equal(RiboCrypt:::mappability_chr_code("chr1", lv), 1L)
  expect_equal(RiboCrypt:::mappability_chr_code("X", lv), 2L)
  expect_equal(RiboCrypt:::mappability_chr_code("chrM", lv), 3L)
  expect_true(is.na(RiboCrypt:::mappability_chr_code("chr9", lv)))
})

test_that("store discovery: newest version next to the genome; option overrides; absent -> NULL", {
  root <- tempfile("ref"); dir.create(file.path(root, "mappability"), recursive = TRUE)
  expect_null(RiboCrypt:::mappability_find_store(file.path(root, "mappability")))
  expect_null(RiboCrypt:::mappability_find_store(file.path(root, "nope")))
  fx <- store_fixture()
  for (v in c("rc_mappability_v2", "rc_mappability_v10")) {
    dir.create(file.path(root, "mappability", v))
    file.copy(list.files(fx, full.names = TRUE), file.path(root, "mappability", v))
  }
  dir.create(file.path(root, "mappability", "rc_mappability_v11"))   # incomplete: ignored
  expect_equal(basename(RiboCrypt:::mappability_find_store(file.path(root, "mappability"))), "rc_mappability_v10")
  withr::local_options(RiboCrypt.mappability_store = c("Mus musculus" = "/nonexistent"))
  expect_null(RiboCrypt:::mappability_store_dir(NULL))
})

test_that("browser uses the coloured strip when a store exists, else the bigWig track", {
  dir <- store_fixture()
  withr::local_options(RiboCrypt.mappability_store = dir)
  display <- GRangesList(GRanges("1", IRanges(c(1001, 1401), c(1100, 1500)), "+"))
  res <- RiboCrypt:::custom_seq_track_panels(NULL, display, phyloP = FALSE, mapability = TRUE)
  expect_named(res$custom_bigwig_panels, "mappability_strip")
  plots <- RiboCrypt:::custom_panels_to_plotly(res$custom_bigwig_panels)
  # strip keeps its hover (not reset by the custom-track axis polish)
  expect_equal(plots$mappability_strip$x$data[[1]]$hoverinfo, "text")
  expect_length(RiboCrypt:::custom_seq_track_panels(NULL, display, FALSE, FALSE)$custom_bigwig_panels, 0)
  # store present but locus outside its domain: no panel (no bigWig fallback)
  far <- GRangesList(GRanges("2", IRanges(1, 100), "+"))
  expect_length(RiboCrypt:::custom_seq_track_panels(NULL, far, FALSE, TRUE)$custom_bigwig_panels, 0)
})

test_that("legend lists every class colour", {
  html <- as.character(RiboCrypt:::mappability_legend_ui())
  for (col in RiboCrypt:::mappability_classes()) expect_match(html, col, fixed = TRUE)
})

test_that("reads select rows by track, part, chromosome and bin (index columns are not shadowed)", {
  dir <- tempfile("rc_map_bins"); dir.create(dir)
  # exact nu: 2 rows in bin 0 of chr 1, 1 row in bin 2 of chr 1, 1 row on chr 2; le2mm nu: 5 rows in bin 0
  ex <- data.table(chr = c(1L, 1L, 1L, 2L), start = c(10L, 50L, 2000010L, 10L), end = c(20L, 60L, 2000020L, 20L),
                   strand = 1L, n = NA_integer_)
  le <- data.table(chr = 1L, start = seq(100L, 500L, 100L), end = seq(110L, 510L, 100L), strand = 0L, n = 2L)
  fst::write_fst(ex, file.path(dir, "exact_nu.fst")); fst::write_fst(le, file.path(dir, "le2mm_nu.fst"))
  fst::write_fst(data.table(track = c("exact", "exact", "exact", "le2mm"), part = "nu", chr = c(1L, 1L, 2L, 1L),
                            bin = c(0L, 2L, 0L, 0L), from = c(1L, 3L, 4L, 1L), to = c(2L, 3L, 4L, 5L)),
                 file.path(dir, "index.fst"))
  saveRDS(list(seqlevels = c("1", "2"), bin = 1e6, tracks = list()), file.path(dir, "meta.rds"))
  store <- RiboCrypt:::mappability_store_open(dir)
  expect_equal(nrow(RiboCrypt:::mappability_read(store, "exact", "nu", 1L, 1L, 1000L)), 2)
  expect_equal(RiboCrypt:::mappability_read(store, "exact", "nu", 1L, 1999990L, 2000030L)$start, 2000010L)
  expect_equal(nrow(RiboCrypt:::mappability_read(store, "exact", "nu", 2L, 1L, 1000L)), 1)
  expect_equal(nrow(RiboCrypt:::mappability_read(store, "le2mm", "nu", 1L, 1L, 1000L)), 5)
  expect_equal(nrow(RiboCrypt:::mappability_read(store, "exact", "nu", 1L, 1L, 2100000L)), 3)
  expect_null(RiboCrypt:::mappability_read(store, "le2mm", "def", 1L, 1L, 1000L))
})
