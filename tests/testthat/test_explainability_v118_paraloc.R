library(GenomicRanges)
library(data.table)

genomic_display2 <- function(tx) GRanges(as.character(seqnames(tx))[1], IRanges(min(start(tx)), max(end(tx))), strand(tx)[1])

# Store of one 200-nt transcript with v1.18 extras: D3 status, dORF type colour, E-term bit mask.
write_v118_store <- function(dir) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  fst::write_fst(data.frame(transcript = "TX1", gene = "G1", n = 200L, cs = 11L, ce = 100L,
                            w_from = 1L, w_to = 3L, l_from = 1L, l_to = 0L, p_from = 1L, p_to = 0L,
                            f_from = 1L, f_to = 0L), file.path(dir, "index.fst"))
  fst::write_fst(data.frame(start = c(1L, 101L, 151L), end = c(100L, 150L, 200L), region = c(2L, 3L, 3L),
                            status = c(1L, 3L, 1L), type = c(NA, 1L, 2L), O = c(10, 30, 50), E = c(10, 12, 48),
                            usable_frac = 1, status_all = 1L, type_all = NA_integer_, O_all = 10, E_all = 10,
                            note = NA_integer_, terms = c(8L, 3L, 5L)), file.path(dir, "windows.fst"))
  fst::write_fst(data.frame(start = integer(), end = integer()), file.path(dir, "lowmap.fst"))
  fst::write_fst(data.frame(codon = integer(), plddt = integer()), file.path(dir, "plddt.fst"))
  fst::write_fst(data.frame(s = integer(), e = integer(), ftype = integer(), note = character()),
                 file.path(dir, "features.fst"))
  saveRDS(list(organism = "Homo sapiens",
               status_levels = c("understood", "RGST", "D3", "not understood"),
               status_colors = c(understood = "#a1d99b", RGST = "#6baed6", D3 = "#dfc27d", "not understood" = "#de2d26"),
               status_labels = c(D3 = "D3 (3'UTR stop-distance profile, descriptive)"),
               type_colors = c("translated ORF (classifier v3 dORF/doORF)" = "#238b45"),
               type_labels = c("translated ORF (classifier v3 dORF/doORF)" = "understood: called dORF/doORF term"),
               term_levels = c("3'UTR RNA term", "D3 stop-distance profile", "called dORF/doORF term",
                               "per-gene CDS contrast + 5' ramp", "Sec UGA step"),
               region_levels = c("leader", "CDS", "trailer"),
               type_levels = c("3'UTR stop-distance profile (D3)", "translated ORF (classifier v3 dORF/doORF)"),
               ftype_levels = "domain", feature_lanes = list("Domains" = "domain")), file.path(dir, "meta.rds"))
  dir
}

write_paraloc_store <- function(dir) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  # TX1: CDS signal peptide 11-79 (UniProt span), overlapping translon compartment 60-130 crossing the junction
  # at tx nt 100 (faint, second sub-lane); TX3 (not in the explain store): MTS of a uORF.
  fst::write_fst(data.frame(transcript = c("TX1", "TX3"), n = c(200L, 50L), cs = c(11L, NA), ce = c(100L, NA),
                            s_from = c(1L, 3L), s_to = c(2L, 3L), gene = c("G1", "G3")), file.path(dir, "index.fst"))
  fst::write_fst(data.frame(start = c(11L, 60L, 5L), end = c(79L, 130L, 40L), kind = c(1L, 6L, 2L),
                            cls = c(3L, 2L, 1L), score = c(0.97, 0.66, 0.81), delta = NA_real_, span = c(1L, 5L, 2L),
                            ttype = c(1L, 3L, 3L), id = c("TX1", "pep7", "pep9"),
                            tools = c("SignalP6 0.99, TargetP2 0.95", "DeepLoc-family localisation winner (loc_prediction)",
                                      "TargetP2 0.80, MitoFates 0.75, MitoProt 0.88"),
                            genomic = c("1:1410-1478:+", "1:1059-1100:+;1:1401-1430:+", "1:2005-2040:+")),
                 file.path(dir, "segments.fst"))
  saveRDS(list(version = "paraloc_track_v1", kind_levels = c("cds_signal", "orf_signal", "nte_gain", "nte_masks",
                                                              "nte_context", "orf_whole"),
               cls_levels = c("MTS", "nucleus", "SP"), span_levels = c("uniprot", "default_len", "extension",
                                                                      "extension+default_len", "whole_orf"),
               ttype_levels = c("CDS", "NTE", "uORF"),
               kind_alpha = c(cds_signal = 1, orf_signal = 1, nte_gain = 1, nte_masks = 0.8, nte_context = 0.35, orf_whole = 0.35),
               cls_colors = c(SP = "#3182bd", MTS = "#e6550d", nucleus = "#756bb1")), file.path(dir, "meta.rds"))
  dir
}

tx1_minus <- GRanges("1", IRanges(c(1401, 1001), c(1500, 1100)), "-")   # tx nt 1-100 = exon 1401-1500

test_that("v1.18 store: D3 status colour, dORF type colour, E-term hover, legend", {
  dir <- write_v118_store(tempfile("explain_v118"))
  store <- RiboCrypt:::explain_store_open(dir)
  explain <- list(active = TRUE, ids = "TX1", layers = "status", tx_exons = GRangesList(TX1 = tx1_minus))
  display <- GRangesList(TX1 = genomic_display2(tx1_minus))
  ex <- RiboCrypt:::explain_locus_data(explain, display, store)
  p <- RiboCrypt:::explain_status_panel(ex, display)
  nm <- vapply(p$x$data, `[[`, "", "name")
  expect_setequal(nm, c("understood", "D3", "translated ORF (classifier v3 dORF/doORF)"))
  col <- setNames(vapply(p$x$data, function(t) t$marker$color, ""), nm)
  expect_equal(col[["D3"]], "#dfc27d")
  expect_equal(col[["translated ORF (classifier v3 dORF/doORF)"]], "#238b45")
  hov <- setNames(lapply(p$x$data, `[[`, "hovertext"), nm)
  expect_match(hov$D3, "<b>D3</b>")
  expect_match(hov$D3, "E terms: 3'UTR RNA term, D3 stop-distance profile")
  expect_match(hov$`translated ORF (classifier v3 dORF/doORF)`, "E terms: 3'UTR RNA term, called dORF/doORF term")
  expect_match(hov$understood, "E terms: per-gene CDS contrast")
  html <- as.character(RiboCrypt:::explain_status_legend_ui(store$meta))
  expect_match(html, "D3 (3'UTR stop-distance profile, descriptive)", fixed = TRUE)
  expect_match(html, "understood: called dORF/doORF term", fixed = TRUE)
  expect_match(html, "#238b45", fixed = TRUE)
  # merged legend over stores keeps the type colours
  html2 <- as.character(RiboCrypt:::explain_status_legend_store(list(store, store)))
  expect_match(html2, "#238b45", fixed = TRUE)
})

test_that("E-term decoding: bits, NA and missing levels", {
  lv <- c("a", "b", "c")
  expect_equal(RiboCrypt:::explain_decode_terms(c(0L, 1L, 6L, NA, 7L), lv), c("", "a", "b, c", "", "a, b, c"))
  expect_equal(RiboCrypt:::explain_decode_terms(c(1L, 2L), NULL), c("", ""))
  expect_equal(RiboCrypt:::explain_decode_terms(NULL, lv), character())
})

test_that("older stores (no terms / type colours) render unchanged", {
  dir <- write_v118_store(tempfile("explain_v117like"))
  m <- readRDS(file.path(dir, "meta.rds")); m$type_colors <- NULL; m$term_levels <- NULL; saveRDS(m, file.path(dir, "meta.rds"))
  store <- RiboCrypt:::explain_store_open(dir)
  explain <- list(active = TRUE, ids = "TX1", layers = "status", tx_exons = GRangesList(TX1 = tx1_minus))
  display <- GRangesList(TX1 = genomic_display2(tx1_minus))
  p <- RiboCrypt:::explain_status_panel(RiboCrypt:::explain_locus_data(explain, display, store), display)
  expect_setequal(vapply(p$x$data, `[[`, "", "name"), c("understood", "D3"))
  expect_false(any(grepl("E terms", unlist(lapply(p$x$data, `[[`, "hovertext")))))
})

test_that("Paralloc lane: discovery next to the explain store, packing, junction split, hover", {
  root <- tempfile("ref_org")
  edir <- write_v118_store(file.path(root, "explainability", "v118"))
  write_paraloc_store(file.path(root, "predicted_translons", "localization", "paraloc_track_v1"))
  dir.create(file.path(root, "predicted_translons", "localization", "paraloc_track_v2"))  # incomplete: ignored
  expect_equal(basename(RiboCrypt:::paraloc_store_dir(edir, "Homo sapiens")), "paraloc_track_v1")
  store <- RiboCrypt:::explain_store_open(edir)
  expect_false(is.null(store$paraloc))
  tl <- RiboCrypt:::explain_stores_translons(list(store))
  expect_true("paraloc" %in% RiboCrypt:::explain_layer_choices(tl))
  # Paralloc-only transcript is a candidate only with the layer selected
  expect_setequal(RiboCrypt:::explain_known_ids(store, c("status", "paraloc")), c("TX1", "TX3"))
  expect_setequal(RiboCrypt:::explain_known_ids(store, "status"), "TX1")
  seg <- RiboCrypt:::paraloc_tx(store$paraloc, "TX1")$seg
  expect_equal(seg$slot, c(1L, 2L))
  explain <- list(active = TRUE, ids = "TX1", layers = c("status", "paraloc"), tx_exons = GRangesList(TX1 = tx1_minus))
  display <- GRangesList(TX1 = genomic_display2(tx1_minus))
  ex <- RiboCrypt:::explain_locus_data(explain, display, store)
  pl <- ex$paraloc
  # compartment segment 60-130 crosses the junction after tx nt 100 -> two boxes
  expect_equal(pl[id == "pep7"]$x0, c(60, 401)); expect_equal(pl[id == "pep7"]$x1, c(100, 430))
  panels <- RiboCrypt:::explain_bottom_panels(ex, display)
  expect_true("explain_paraloc" %in% names(panels))
  tr <- panels$explain_paraloc$x$data
  nm <- vapply(tr, `[[`, "", "name")
  expect_setequal(nm, c("SP", "nucleus"))
  sp <- tr[[which(nm == "SP")]]
  expect_match(sp$hovertext, "signal peptide \\(SP\\)")
  expect_match(sp$hovertext, "predictor: SignalP6 \\+ TargetP2 \\(mean\\): SignalP6 0.99, TargetP2 0.95")
  expect_match(sp$hovertext, "score 0.97")
  expect_match(sp$hovertext, "UniProt span")
  expect_equal(as.numeric(sp$marker$opacity), 1)
  nuc <- tr[[which(nm == "nucleus")]]
  expect_equal(as.numeric(nuc$marker$opacity), c(0.35, 0.35))
  expect_match(nuc$hovertext[1], "predicted compartment")
  expect_gt(as.numeric(sp$base), as.numeric(nuc$base[1]))
  expect_equal(panels$explain_paraloc$x$layout$yaxis$range, c(0, 2))
  html <- as.character(RiboCrypt:::paraloc_legend_ui(store$paraloc$meta))
  expect_match(html, "paraloc_track_v1"); expect_match(html, "#e6550d")
  # length mismatch: skipped with a warning
  bad <- GRanges("1", IRanges(2001, 2020), "+")
  expect_warning(res <- RiboCrypt:::paraloc_locus(store$paraloc, "TX3", bad, bad), "length mismatch")
  expect_null(res)
})

test_that("Paralloc store option overrides discovery; no store -> no layer", {
  dir <- write_v118_store(tempfile("explain_nopl"))
  store <- RiboCrypt:::explain_store_open(dir)
  expect_null(store$paraloc)
  expect_false("paraloc" %in% RiboCrypt:::explain_layer_choices(RiboCrypt:::explain_stores_translons(list(store))))
  pdir <- write_paraloc_store(tempfile("pl"))
  withr::local_options(RiboCrypt.paraloc_store = c("Homo sapiens" = pdir))
  expect_equal(RiboCrypt:::paraloc_store_dir(dir, "Homo sapiens"), normalizePath(pdir))
  expect_null(RiboCrypt:::paraloc_store_dir(dir, "Mus musculus"))
})

test_that("interval packing never overlaps within a lane", {
  s <- c(1L, 5L, 12L, 3L, 20L); e <- c(10L, 15L, 18L, 4L, 25L)
  sl <- RiboCrypt:::paraloc_pack(s, e)
  for (k in unique(sl)) {
    i <- which(sl == k); o <- order(s[i])
    expect_true(all(s[i][o][-1] > e[i][o][-length(i)]))
  }
  expect_length(RiboCrypt:::paraloc_pack(integer(), integer()), 0)
})
