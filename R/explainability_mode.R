# Explainability mode for the single-locus browser.
#
# Shows precomputed explainer layers for the transcripts of an organism store (human MANE, mouse canonical):
#   - window status ribbon (90-nt windows, status_fam, hover: type, O, E, log2 O/E)
#   - expected P-site coverage E overlaid on Ribo-seq tracks (scaled to each library)
#   - AlphaFold pLDDT per codon
#   - UniProt features per protein (domains, disorder, topology, sites, AF helix/strand)
#   - translon calls (optional, private classifier release; any transcript of the release, not only MANE)
# Data come from a compact store directory (built at release time by
# build_explain_store.R): index.fst (one row per transcript with row ranges),
# windows.fst, lowmap.fst, plddt.fst, features.fst and meta.rds. Rows are read
# lazily per transcript and kept in an in-memory cache. Optional extra lane files
# (translons_<version>.fst, translons_<version>_index.fst, translons_<version>_meta.rds)
# add the translon calls lane; stores without them work as before.
#
# Stores are given by option "RiboCrypt.explain_store" or the environment variable
# RIBOCRYPT_EXPLAIN_STORE; without them the mode is hidden. Accepted values:
#   - one store directory (organism = meta.rds$organism), as before
#   - a named character vector / list, names = organism, e.g.
#     c("Homo sapiens" = "<dir>/v117r", "Mus musculus" = "<dir>/v117m")
#   - an unnamed vector of store directories (organism of each from its meta.rds)
#   - environment variable: one directory, or "Homo sapiens=<dir>;Mus musculus=<dir>"
# The store of the selected experiment's organism is used.

#' @noRd
explain_layer_choices <- function(store = NULL) {
  ch <- c("Window status" = "status", "Expected E on P-sites" = "expected",
          "AlphaFold pLDDT" = "plddt", "UniProt features" = "features")
  if (!is.null(store$translons)) ch <- c(ch, "Translon calls" = "translons")
  if (isTRUE(store$has_paraloc)) ch <- c(ch, "Localisation signals (Paralloc)" = "paraloc")
  c(ch, "Exons (alternating)" = "exons")
}

#' Normalise an organism name for matching ("Homo_sapiens" == "homo sapiens").
#' @noRd
explain_norm_organism <- function(x) tolower(gsub("[_ ]+", " ", trimws(as.character(x))))

#' Raw store configuration: option, else environment variable (parsed).
#' @noRd
explain_store_config <- function() {
  d <- getOption("RiboCrypt.explain_store", NULL)
  if (is.null(d)) {
    env <- Sys.getenv("RIBOCRYPT_EXPLAIN_STORE", "")
    if (!nzchar(env)) return(NULL)
    parts <- trimws(strsplit(env, ";", fixed = TRUE)[[1]])
    parts <- parts[nzchar(parts)]
    named <- grepl("=", parts, fixed = TRUE)
    d <- ifelse(named, sub("^[^=]*=", "", parts), parts)
    names(d) <- ifelse(named, trimws(sub("=.*$", "", parts)), "")
    if (!any(named)) d <- unname(d)
  }
  d
}

#' Configured explainability stores: named character vector (names = normalised
#' organism, values = store directory), empty when none is configured / valid.
#'
#' Unnamed entries take the organism from the store's meta.rds. Directories
#' without index.fst are dropped; for duplicated organisms the first entry wins.
#' @noRd
explain_store_dirs <- function(config = explain_store_config()) {
  empty <- setNames(character(), character())
  if (is.null(config) || !length(config)) return(empty)
  nms <- names(config)
  if (is.null(nms)) nms <- rep("", length(config))
  dirs <- vapply(config, function(x) if (length(x) == 1 && !is.na(x)) as.character(x) else "", "")
  ok <- nzchar(dirs) & dir.exists(dirs) & file.exists(file.path(dirs, "index.fst"))
  if (!any(ok)) return(empty)
  dirs <- normalizePath(dirs[ok]); nms <- nms[ok]
  orgs <- vapply(seq_along(dirs), function(i) {
    if (!is.na(nms[i]) && nzchar(nms[i])) return(explain_norm_organism(nms[i]))
    explain_store_organism(dirs[i])
  }, "")
  keep <- !is.na(orgs) & nzchar(orgs) & !duplicated(orgs)
  setNames(unname(dirs[keep]), orgs[keep])
}

#' Organism (normalised) of a store directory from its meta.rds (cached), "" when unreadable.
#' @noRd
explain_store_organism <- function(dir) {
  key <- paste0("org:", dir)
  org <- .rc_explain_env[[key]]
  if (!is.null(org)) return(org)
  org <- tryCatch(explain_norm_organism(readRDS(file.path(dir, "meta.rds"))$organism),
                  error = function(e) "")
  if (length(org) != 1 || is.na(org)) org <- ""
  assign(key, org, envir = .rc_explain_env)
  org
}

#' Explainability store directory, NULL when not configured.
#'
#' @param organism organism name of the selected experiment; NULL = the first
#'  configured store (old single-store behaviour, and "is any store configured").
#' @noRd
explain_store_dir <- function(organism = NULL) {
  dirs <- explain_store_dirs()
  if (!length(dirs)) return(NULL)
  if (is.null(organism)) return(unname(dirs[1]))
  if (length(organism) != 1 || is.na(organism)) return(NULL)
  d <- dirs[explain_norm_organism(organism)]
  if (is.na(d)) NULL else unname(d)
}

#' Open the store of an organism, NULL when that organism has none.
#' @noRd
explain_store_for <- function(organism) {
  dir <- explain_store_dir(organism)
  if (is.null(dir)) return(NULL)
  explain_store_open(dir)
}

#' All configured stores, opened (unreadable ones dropped).
#' @noRd
explain_stores_open <- function() {
  dirs <- explain_store_dirs()
  stores <- lapply(dirs, function(d) tryCatch(explain_store_open(d), error = function(e) NULL))
  stores[!vapply(stores, is.null, logical(1))]
}

#' Short display names of the configured organisms ("human, mouse").
#' @noRd
explain_organism_labels <- function(orgs = names(explain_store_dirs())) {
  short <- c("homo sapiens" = "human", "mus musculus" = "mouse")
  lab <- ifelse(orgs %in% names(short), short[orgs], orgs)
  paste(unname(lab), collapse = ", ")
}

.rc_explain_env <- new.env(parent = emptyenv())

#' Open (and cache) an explainability store.
#' @noRd
explain_store_open <- function(dir = explain_store_dir()) {
  if (is.null(dir)) return(NULL)
  key <- paste0("store:", dir)
  store <- .rc_explain_env[[key]]
  if (!is.null(store)) return(store)
  index <- fst::read_fst(file.path(dir, "index.fst"), as.data.table = TRUE)
  meta <- readRDS(file.path(dir, "meta.rds"))
  store <- list(dir = dir, index = index, meta = meta,
                row = setNames(seq_len(nrow(index)), index$transcript),
                cache = new.env(parent = emptyenv()),
                translons = explain_translons_open(dir),
                paraloc = paraloc_open(paraloc_store_dir(dir, meta$organism)))
  assign(key, store, envir = .rc_explain_env)
  store
}

#' Translon calls lane files of a store (newest version), NULL when absent.
#'
#' Files: translons_<v>.fst (one row per call, rows contiguous per transcript),
#' translons_<v>_index.fst (transcript, gene, n, chr, strand, from, to, n_slots)
#' and translons_<v>_meta.rds (levels, colours). Only the index and meta are read
#' here; call rows are read per transcript (explain_translons_tx).
#' @noRd
explain_translons_open <- function(dir) {
  idx <- list.files(dir, pattern = "^translons_v[0-9_]+_index\\.fst$")
  if (!length(idx)) return(NULL)
  vers <- sub("^translons_(v[0-9_]+)_index\\.fst$", "\\1", idx)
  nv <- numeric_version(gsub("_", ".", sub("^v", "", vers)), strict = FALSE)
  v <- vers[order(nv, decreasing = TRUE)][1]
  files <- file.path(dir, paste0("translons_", v, c(".fst", "_index.fst", "_meta.rds")))
  if (!all(file.exists(files))) return(NULL)
  res <- tryCatch({
    index <- fst::read_fst(files[2], as.data.table = TRUE)
    meta <- readRDS(files[3])
    list(file = files[1], version = v, index = index, meta = meta,
         row = setNames(seq_len(nrow(index)), index$transcript),
         cache = new.env(parent = emptyenv()))
  }, error = function(e) {
    warning("Explainability: could not read translon lane ", v, ": ", conditionMessage(e))
    NULL
  })
  res
}

#' Translon calls of one transcript (lazy, cached); NULL when none.
#' @noRd
explain_translons_tx <- function(tl, id, max_cache = 2000L) {
  if (is.null(tl)) return(NULL)
  hit <- tl$cache[[id]]
  if (!is.null(hit)) return(hit)
  i <- tl$row[id]
  if (is.na(i)) return(NULL)
  r <- tl$index[i]
  calls <- fst::read_fst(tl$file, from = r$from, to = r$to, as.data.table = TRUE)
  res <- list(tx = r, calls = calls)
  if (length(ls(tl$cache)) >= max_cache) rm(list = ls(tl$cache), envir = tl$cache)
  assign(id, res, envir = tl$cache)
  res
}

#' Transcript ids with any data for the requested layers.
#' @noRd
explain_known_ids <- function(store, layers = NULL) {
  ids <- names(store$row)
  if (!is.null(store$translons) && "translons" %in% layers) {
    ids <- union(ids, names(store$translons$row))
  }
  if (!is.null(store$paraloc) && "paraloc" %in% layers) {
    ids <- union(ids, names(store$paraloc$row))
  }
  ids
}

#' Is the organism supported?
#'
#' @param store a store (check that store only); NULL = any configured store.
#' @noRd
explain_organism_supported <- function(organism_name, store = NULL) {
  if (length(organism_name) != 1 || is.na(organism_name)) return(FALSE)
  if (is.null(store)) return(!is.null(explain_store_dir(organism_name)))
  identical(explain_norm_organism(organism_name), explain_norm_organism(store$meta$organism))
}

#' Read all layers of one transcript (lazy, cached).
#' @noRd
explain_store_tx <- function(store, id, max_cache = 2000L) {
  hit <- store$cache[[id]]
  if (!is.null(hit)) return(hit)
  i <- store$row[id]
  if (is.na(i)) return(NULL)
  r <- store$index[i]
  rd <- function(file, from, to) {
    if (to < from) return(NULL)
    fst::read_fst(file.path(store$dir, file), from = from, to = to, as.data.table = TRUE)
  }
  res <- list(tx = r,
              windows = rd("windows.fst", r$w_from, r$w_to),
              lowmap = rd("lowmap.fst", r$l_from, r$l_to),
              plddt = rd("plddt.fst", r$p_from, r$p_to),
              features = rd("features.fst", r$f_from, r$f_to))
  if (length(ls(store$cache)) >= max_cache) rm(list = ls(store$cache), envir = store$cache)
  assign(id, res, envir = store$cache)
  res
}

#' Transcript ids of the store to show for the current locus.
#'
#' The selected transcript, plus (genomic view with full annotation) every store
#' transcript overlapping the displayed region.
#' @noRd
explain_candidate_ids <- function(store, selected, tx_grl, display_region,
                                  include_overlapping = FALSE, max_tx = 5L, layers = NULL) {
  if (is.null(store)) return(character())
  known <- explain_known_ids(store, layers)
  ids <- intersect(selected, known)
  if (include_overlapping) {
    spans <- explain_store_spans(store, tx_grl, known)
    if (length(spans)) {
      hits <- names(spans)[IRanges::overlapsAny(spans, unlist(display_region, use.names = FALSE),
                                       ignore.strand = TRUE)]
      ids <- unique(c(ids, hits))
    }
  }
  ids <- intersect(ids, names(tx_grl))
  head(ids, max_tx)
}

#' Genomic spans of store transcripts in one annotation (cached per annotation).
#' @noRd
explain_store_spans <- function(store, tx_grl, known = names(store$row)) {
  key <- paste0("spans:", length(tx_grl), ":", names(tx_grl)[1], ":",
                names(tx_grl)[length(tx_grl)], ":", length(known))
  spans <- store$cache[[key]]
  if (is.null(spans)) {
    ids <- intersect(known, names(tx_grl))
    spans <- unlist(BiocGenerics::range(tx_grl[ids]))
    assign(key, spans, envir = store$cache)
  }
  spans
}

# ---- coordinate mapping ---------------------------------------------------

#' Linear segments mapping transcript nt to browser display x.
#'
#' The browser x axis is the 1-based position along the tiled display region
#' (5' to 3' in display direction: the transcript in transcript view, the gene
#' span or pseudo-intron exons in genomic view). Each exon part inside the
#' display becomes one segment: transcript nt tx_a..tx_b maps to x_a..x_b, with
#' x_b = x_a + sgn * (tx_b - tx_a).
#' @param tx_exons GRanges, exons of one transcript in transcript order (ORFik convention)
#' @param display GRanges, display region in display order
#' @return data.table(tx_a, tx_b, x_a, sgn), ordered by tx_a
#' @noRd
explain_tx_display_map <- function(tx_exons, display) {
  empty <- data.table(tx_a = integer(), tx_b = integer(), x_a = integer(), sgn = integer())
  if (!length(tx_exons) || !length(display)) return(empty)
  tx_minus <- as.character(strand(tx_exons))[1] == "-"
  ts0 <- start(tx_exons); te0 <- end(tx_exons)
  ew <- te0 - ts0 + 1L
  etx0 <- cumsum(c(0L, ew[-length(ew)])) + 1L
  ds0 <- start(display); de0 <- end(display)
  dw <- de0 - ds0 + 1L
  dcum0 <- cumsum(c(0L, dw[-length(dw)]))
  # both sets are non-overlapping: intersect on genomic coordinates sorted by start
  oq <- order(ts0); od <- order(ds0)
  hit <- explain_intersect_sorted(ts0[oq], te0[oq], ds0[od], de0[od])
  if (!length(hit$qi)) return(empty)
  q <- oq[hit$qi]; s <- od[hit$si]
  ps <- hit$lo; pe <- hit$hi
  if (tx_minus) {
    tx_a <- etx0[q] + (te0[q] - pe); g_a <- pe
  } else {
    tx_a <- etx0[q] + (ps - ts0[q]); g_a <- ps
  }
  tx_b <- tx_a + (pe - ps)
  disp_minus <- as.character(strand(display))[1] == "-"
  x_a <- if (disp_minus) dcum0[s] + de0[s] - g_a + 1L else dcum0[s] + g_a - ds0[s] + 1L
  sgn <- if (tx_minus == disp_minus) 1L else -1L
  o <- order(tx_a)
  data.table(tx_a = as.integer(tx_a[o]), tx_b = as.integer(tx_b[o]), x_a = as.integer(x_a[o]),
             sgn = rep(as.integer(sgn), length(o)))
}

#' Intersect query ranges with sorted, non-overlapping segments.
#'
#' Plain vector arithmetic (findInterval), no S4 / overlap index overhead.
#' @return list(qi, si, lo, hi): query index, segment index, intersection bounds
#' @noRd
explain_intersect_sorted <- function(a, b, sa, sb) {
  k1 <- findInterval(a - 1L, sb) + 1L  # first segment ending at or after a
  k2 <- findInterval(b, sa)            # last segment starting at or before b
  n <- pmax(k2 - k1 + 1L, 0L)
  keep <- n > 0L
  qi <- rep.int(seq_along(a), n)
  si <- sequence(n[keep], from = k1[keep])
  list(qi = qi, si = si, lo = pmax(a[qi], sa[si]), hi = pmin(b[qi], sb[si]))
}

#' Map transcript ranges to display x ranges, split at exon junctions.
#'
#' Pieces of one range that stay contiguous on the display (transcript view)
#' are merged again.
#' @return data.table(id, x0, x1) with x0 <= x1; id = index of the input range
#' @noRd
explain_map_ranges <- function(start, end, map) {
  empty <- data.table(id = integer(), x0 = integer(), x1 = integer())
  if (!length(start) || !nrow(map)) return(empty)
  hit <- explain_intersect_sorted(as.integer(start), as.integer(end), map$tx_a, map$tx_b)
  if (!length(hit$qi)) return(empty)
  si <- hit$si
  xl <- map$x_a[si] + map$sgn[si] * (hit$lo - map$tx_a[si])
  xh <- map$x_a[si] + map$sgn[si] * (hit$hi - map$tx_a[si])
  id <- hit$qi; x0 <- pmin(xl, xh); x1 <- pmax(xl, xh)
  o <- order(id, x0)
  id <- id[o]; x0 <- x0[o]; x1 <- x1[o]
  n <- length(id)
  if (n > 1L) {
    new_grp <- c(TRUE, id[-1L] != id[-n] | x0[-1L] != x1[-n] + 1L)
    g <- cumsum(new_grp)
    first <- which(new_grp); last <- c(first[-1L] - 1L, n)
    return(data.table(id = id[first], x0 = x0[first], x1 = x1[last]))
  }
  data.table(id = id, x0 = x0, x1 = x1)
}

#' Map transcript nt positions to display x (NA when outside the display).
#' @noRd
explain_map_points <- function(pos, map) {
  if (!nrow(map)) return(rep(NA_integer_, length(pos)))
  k <- findInterval(pos, map$tx_a)
  ok <- k > 0L
  ok[ok] <- pos[ok] <= map$tx_b[k[ok]]
  x <- rep(NA_integer_, length(pos))
  x[ok] <- map$x_a[k[ok]] + map$sgn[k[ok]] * (pos[ok] - map$tx_a[k[ok]])
  x
}

# ---- layer data -----------------------------------------------------------

#' Build all requested layers for the current locus in display coordinates.
#' @param explain list from the controller: active, layers, ids, tx_exons
#' @param display_range GRangesList (one group), the final browser display region
#' @noRd
explain_locus_data <- function(explain, display_range,
                               store = explain_store_open(if (is.null(explain$store_dir)) explain_store_dir()
                                                          else explain$store_dir)) {
  if (is.null(store) || !isTRUE(explain$active) || !length(explain$ids)) return(NULL)
  layers <- explain$layers
  display <- unlist(display_range, use.names = FALSE)
  out <- list(status = list(), expected = list(), plddt = list(), features = list(),
              translons = list(), paraloc = list(), meta = store$meta, ids = character())
  for (id in explain$ids) {
    if (!(id %in% names(explain$tx_exons))) next
    exons <- explain$tx_exons[[id]]
    if ("translons" %in% layers && !is.null(store$translons)) {
      out$translons[[id]] <- explain_translons_locus(store$translons, id, exons, display)
    }
    if ("paraloc" %in% layers && !is.null(store$paraloc)) {
      out$paraloc[[id]] <- paraloc_locus(store$paraloc, id, exons, display)
    }
    d <- explain_store_tx(store, id)
    if (is.null(d)) next
    if (sum(width(exons)) != d$tx$n) {
      warning("Explainability: transcript length mismatch for ", id, ", skipped")
      next
    }
    map <- explain_tx_display_map(exons, display)
    if (!nrow(map)) next
    out$ids <- c(out$ids, id)
    slane <- length(out$ids)  # store lanes stay contiguous when ids without store data are skipped
    label <- paste0(d$tx$gene, " (", id, ")")
    if (any(c("status", "expected") %in% layers) && !is.null(d$windows)) {
      w <- d$windows
      if ("status" %in% layers) {
        p <- explain_map_ranges(w$start, w$end, map)
        if (nrow(p)) out$status[[id]] <- data.table(p, w[p$id], lane = slane, label = label)
      }
      if ("expected" %in% layers) out$expected[[id]] <- explain_expected_segments(d, map)
    }
    if ("plddt" %in% layers && !is.null(d$plddt)) {
      out$plddt[[id]] <- explain_plddt_bins(d, map, lane = slane, label = label)
    }
    if ("features" %in% layers && !is.null(d$features)) {
      f <- d$features
      cs <- d$tx$cs
      p <- explain_map_ranges(cs + 3L * (f$s - 1L), cs + 3L * f$e - 1L, map)
      if (nrow(p)) out$features[[id]] <- data.table(p, f[p$id], label = label)
    }
  }
  out$translons <- rbindlist(out$translons)
  if (nrow(out$translons)) out$translons_meta <- store$translons$meta
  out$paraloc <- rbindlist(out$paraloc)
  if (nrow(out$paraloc)) out$paraloc_meta <- store$paraloc$meta
  if (!length(out$ids) && !nrow(out$translons) && !nrow(out$paraloc)) return(NULL)
  out$status <- rbindlist(out$status)
  out$expected <- rbindlist(out$expected)
  out$plddt <- rbindlist(out$plddt)
  out$features <- rbindlist(out$features)
  out
}

#' Translon calls of one transcript in display coordinates.
#'
#' ORFs crossing an exon junction become several boxes (genomic view) with the
#' same hover; slot = precomputed sub-lane (packing in transcript nt).
#' @return data.table(x0, x1, call columns, slot, n_slots, tx, label) or NULL
#' @noRd
explain_translons_locus <- function(tl, id, exons, display) {
  d <- explain_translons_tx(tl, id)
  if (is.null(d) || !nrow(d$calls)) return(NULL)
  if (sum(width(exons)) != d$tx$n) {
    warning("Explainability: transcript length mismatch for ", id, " (translon lane), skipped")
    return(NULL)
  }
  map <- explain_tx_display_map(exons, display)
  if (!nrow(map)) return(NULL)
  p <- explain_map_ranges(d$calls$start, d$calls$end, map)
  if (!nrow(p)) return(NULL)
  data.table(x0 = p$x0, x1 = p$x1, d$calls[p$id], n_slots = d$tx$n_slots, tx = id,
             label = paste0(d$tx$gene, " (", id, ")"))
}

#' Per-nt expected density (unique track) at usable positions, in display x.
#'
#' Density of window w = E_w / (width_w * usable_frac_w), the same expected
#' profile as the explainer plot cache (exp_win / 3), undefined in
#' low-mappability intervals.
#' @return data.table(x0, x1, density, E_window)
#' @noRd
explain_expected_segments <- function(d, map) {
  w <- d$windows
  n <- d$tx$n
  # usable = 1..n minus low-mappability intervals (sorted, merged here)
  if (!is.null(d$lowmap) && nrow(d$lowmap)) {
    ls_ <- d$lowmap$start; le_ <- d$lowmap$end
    o <- order(ls_); ls_ <- ls_[o]; le_ <- le_[o]
    cm <- cummax(le_)
    brk <- c(TRUE, ls_[-1L] > cm[-length(cm)] + 1L)
    g <- cumsum(brk)
    ms <- ls_[brk]; me <- as.integer(tapply(le_, g, max))
    us <- c(1L, me + 1L); ue <- c(ms - 1L, n)
    ok <- ue >= us
    us <- us[ok]; ue <- ue[ok]
  } else {
    us <- 1L; ue <- n
  }
  hit <- explain_intersect_sorted(us, ue, w$start, w$end)
  if (!length(hit$qi)) return(NULL)
  k <- hit$si
  dens <- w$E[k] / ((w$end[k] - w$start[k] + 1) * pmax(w$usable_frac[k], 0.05))
  p <- explain_map_ranges(hit$lo, hit$hi, map)
  if (!nrow(p)) return(NULL)
  shown <- unique(k[p$id])
  data.table(x0 = p$x0, x1 = p$x1, density = dens[p$id], window = k[p$id],
             E_window = w$E[k[p$id]], E_shown_total = sum(w$E[shown]))
}

#' pLDDT binned along the CDS (exact per codon for proteins up to max_bins codons).
#' @noRd
explain_plddt_bins <- function(d, map, lane = 1L, label = "", max_bins = 1500L) {
  cod <- d$plddt$codon
  v <- d$plddt$plddt / 10
  b <- max(1L, as.integer(ceiling(max(cod) / max_bins)))
  if (b > 1L) {
    bin <- (cod - 1L) %/% b
    first <- which(!duplicated(bin)); last <- c(first[-1L] - 1L, length(bin))
    cv <- c(0, cumsum(v))
    c1 <- cod[first]; c2 <- cod[last]; v <- (cv[last + 1L] - cv[first]) / (last - first + 1L)
  } else {
    c1 <- c2 <- cod
  }
  cs <- d$tx$cs
  p <- explain_map_ranges(cs + 3L * (c1 - 1L), cs + 3L * c2 - 1L, map)
  if (!nrow(p)) return(NULL)
  data.table(p, plddt = v[p$id], c1 = c1[p$id], c2 = c2[p$id], lane = lane, label = label)
}

# ---- plotly panels --------------------------------------------------------

#' AlphaFold pLDDT confidence bands and colours.
#' @noRd
explain_plddt_band <- function(x) {
  cut(x, c(-Inf, 50, 70, 90, Inf), right = FALSE,
      labels = c("very low (<50)", "low (50-70)", "confident (70-90)", "very high (>90)"))
}
explain_plddt_colors <- function() {
  c("very low (<50)" = "#FF7D45", "low (50-70)" = "#FFDB13",
    "confident (70-90)" = "#65CBF3", "very high (>90)" = "#0053D6")
}

#' Default window status colours and legend labels.
#'
#' "RGST" = reproducible gene-specific trend: windows that are not understood
#' under the mechanistic model and are explained only by per-gene empirical shape
#' terms (per-gene contrast + 5' ramp; explainer v1.15).
#' @noRd
explain_status_defaults <- function() {
  list(colors = c("understood" = "#a1d99b", "RGST" = "#6baed6", "not understood" = "#de2d26",
                  "study-specific" = "#fd8d3c", "family-level" = "#9e9ac8", "unknown" = "#d9d9d9",
                  "contamination" = "#8c510a"),
       labels = c("RGST" = "RGST (reproducible gene-specific trend)"))
}

#' Status colours for a store: store colours override the defaults; any status
#' level of the store without a colour gets a neutral fallback, so windows with
#' a new status value are still drawn.
#' @noRd
explain_status_colors <- function(meta) {
  cols <- explain_status_defaults()$colors
  sc <- meta$status_colors
  if (length(sc)) cols[names(sc)] <- sc
  miss <- setdiff(meta$status_levels, names(cols))
  if (length(miss)) cols[miss] <- "#bdbdbd"
  lv <- meta$status_levels
  if (length(lv)) cols <- cols[c(intersect(lv, names(cols)), setdiff(names(cols), lv))]
  cols
}

#' Legend labels per status (store meta$status_labels override the defaults).
#' @noRd
explain_status_labels <- function(meta) {
  lv <- if (length(meta$status_levels)) meta$status_levels else names(explain_status_defaults()$colors)
  lab <- setNames(lv, lv)
  def <- explain_status_defaults()$labels
  def <- def[intersect(names(def), lv)]
  lab[names(def)] <- def
  sl <- meta$status_labels
  if (length(sl)) lab[intersect(names(sl), lv)] <- sl[intersect(names(sl), lv)]
  lab
}

#' Colour legend of the window status ribbon (HTML).
#' @noRd
explain_status_legend_ui <- function(meta) {
  cols <- explain_status_colors(meta)
  lab <- explain_status_labels(meta)
  lv <- names(lab)
  items <- lapply(lv, function(s) tags$span(
    style = "display:inline-block; margin-right:10px; white-space:nowrap;",
    tags$span(style = paste0("display:inline-block; width:12px; height:12px; margin-right:4px; ",
                             "vertical-align:middle; background:", cols[[s]], ";")),
    lab[[s]]))
  tc <- meta$type_colors
  if (length(tc)) {
    tl <- meta$type_labels
    items <- c(items, lapply(names(tc), function(s) tags$span(
      style = "display:inline-block; margin-right:10px; white-space:nowrap;",
      tags$span(style = paste0("display:inline-block; width:12px; height:12px; margin-right:4px; ",
                               "vertical-align:middle; background:", tc[[s]], ";")),
      if (!is.null(tl) && !is.na(tl[s])) tl[[s]] else s)))
  }
  tags$div(class = "explain-status-legend", style = "font-size:12px;", tags$b("Window status: "), items)
}

#' Feature lane colours.
#' @noRd
explain_feature_lane_colors <- function() {
  c("Domains" = "#1b9e77", "Disorder / bias" = "#d95f02", "Topology" = "#7570b3",
    "Sites / PTM" = "#e7298a", "Helix / strand (AF)" = "#66a61e")
}

#' One bar trace per category; bars span x0..x1 on the shared browser x axis.
#'
#' Traces are plain lists (already in plotly.js form), so the panel never goes
#' through plotly_build (the slow part for a few thousand bars).
#' @noRd
explain_bar_traces <- function(dt, category, colors, y, base, hover) {
  traces <- list()
  for (cat_name in intersect(names(colors), unique(category))) {
    sel <- which(category == cat_name)
    traces[[length(traces) + 1L]] <- list(
      type = "bar", x = I((dt$x0[sel] + dt$x1[sel]) / 2), width = I(dt$x1[sel] - dt$x0[sel] + 1),
      y = I(y[sel]), base = I(base[sel]), name = cat_name, hovertext = I(hover[sel]),
      hoverinfo = "text", textposition = "none", showlegend = FALSE, orientation = "v",
      marker = list(color = unname(colors[cat_name]), line = list(width = 0)),
      xaxis = "x", yaxis = "y")
  }
  traces
}

#' Empty, already built plotly object used as panel template (cached).
#' @noRd
explain_plotly_template <- function() {
  tpl <- .rc_explain_env$template
  if (is.null(tpl)) {
    tpl <- suppressWarnings(plotly::plotly_build(plotly::plot_ly(type = "bar", x = 1, y = 1)))
    tpl$x$data <- list()
    tpl$x$attrs <- list()
    tpl$x$layoutAttrs <- list()
    assign("template", tpl, envir = .rc_explain_env)
  }
  tpl
}

#' Materialized explainability panel on the browser x axis.
#' @noRd
explain_panel <- function(traces, display_range, ylab, yrange) {
  p <- explain_plotly_template()
  # Explain panels are short (~50 px): a rotated y title (22 px) is longer than the
  # panel and runs into the title of the neighbouring panel. Use a small horizontal
  # label in the left margin instead (annotation in the panel's y domain).
  yaxis <- custom_seq_track_yaxis("")
  yaxis$range <- yrange
  yaxis$autorange <- FALSE
  xaxis <- custom_seq_track_xaxis(display_range)
  xaxis$autorange <- FALSE
  p$x$data <- traces
  p$x$layout$xaxis <- xaxis
  p$x$layout$yaxis <- yaxis
  p$x$layout$annotations <- list(explain_panel_label(ylab))
  p$x$layout$bargap <- 0
  p$x$layout$barmode <- "overlay"
  p$x$layout$margin <- list(t = 0, r = 0, b = 0, l = 0, pad = 0)
  p$x$layout$showlegend <- FALSE
  p
}

#' Horizontal panel label left of the plot area, vertically centred on the panel.
#' @noRd
explain_panel_label <- function(text) {
  list(text = text, xref = "paper", x = 0, xanchor = "right", xshift = -4,
       yref = "y domain", y = 0.5, yanchor = "middle", showarrow = FALSE,
       font = list(size = 11))
}

#' Window status ribbon (one lane per transcript).
#' @noRd
explain_status_panel <- function(ex, display_range) {
  dt <- ex$status
  meta <- ex$meta
  status <- meta$status_levels[dt$status]
  type <- meta$type_levels[dt$type]
  status_all <- meta$status_levels[dt$status_all]
  region <- meta$region_levels[dt$region]
  # optional hover-only notes (store v1.17+): integer column `note` decoded by meta$note_levels
  note <- if (!is.null(dt$note) && !is.null(meta$note_levels)) meta$note_levels[dt$note] else rep(NA_character_, nrow(dt))
  hover <- sprintf(
    "%s<br>%s window %d-%d nt: <b>%s</b>%s%s<br>O = %s, E = %s, log2 O/E = %.2f<br>all-mapper track: %s",
    dt$label, region, dt$start, dt$end, status,
    ifelse(is.na(type), "", paste0("<br>", type)),
    ifelse(is.na(note), "", paste0("<br><i>", note, "</i>")),
    format(round(dt$O), big.mark = ",", trim = TRUE), format(round(dt$E), big.mark = ",", trim = TRUE),
    log2((dt$O + 1) / (dt$E + 1)), ifelse(is.na(status_all), "NA", status_all))
  # optional E-term bit mask (store v1.18+): integer column `terms`, levels meta$term_levels (bit k = 2^(k-1))
  terms <- explain_decode_terms(dt$terms, meta$term_levels)
  hover <- paste0(hover, ifelse(nzchar(terms), paste0("<br>E terms: ", terms), ""))
  # optional per-type ribbon colours (store v1.18+: meta$type_colors, e.g. the called dORF term)
  category <- status
  tc <- meta$type_colors
  if (length(tc)) {
    hit <- !is.na(type) & type %in% names(tc)
    category[hit] <- type[hit]
  }
  cols <- c(explain_status_colors(meta), tc[setdiff(names(tc), names(explain_status_colors(meta)))])
  nl <- max(dt$lane)
  traces <- explain_bar_traces(dt, category, cols, rep(0.9, nrow(dt)),
                               nl - dt$lane + 0.05, hover)
  explain_panel(traces, display_range, "Explain", c(0, nl))
}

#' Decode the E-term bit mask of windows ("3'UTR RNA term, D3 ..."); "" when none / no levels.
#' @noRd
explain_decode_terms <- function(terms, levels) {
  n <- length(terms)
  if (is.null(terms) || !length(levels) || !n) return(rep("", n))
  out <- rep("", n)
  terms[is.na(terms)] <- 0L
  for (k in seq_along(levels)) {
    on <- bitwAnd(as.integer(terms), as.integer(2^(k - 1))) > 0L
    out[on] <- ifelse(nzchar(out[on]), paste0(out[on], ", ", levels[k]), levels[k])
  }
  out
}

#' AlphaFold pLDDT panel (bar height = pLDDT, colour = confidence band).
#' @noRd
explain_plddt_panel <- function(ex, display_range) {
  dt <- ex$plddt
  band <- as.character(explain_plddt_band(dt$plddt))
  aa <- ifelse(dt$c1 == dt$c2, paste0("residue ", dt$c1), paste0("residues ", dt$c1, "-", dt$c2))
  hover <- sprintf("%s<br>%s: pLDDT %.1f (%s)", dt$label, aa, dt$plddt, band)
  traces <- explain_bar_traces(dt, band, explain_plddt_colors(), dt$plddt, rep(0, nrow(dt)), hover)
  explain_panel(traces, display_range, "pLDDT", c(0, 100))
}

#' UniProt feature panel (one lane per feature group).
#' @noRd
explain_features_panel <- function(ex, display_range) {
  dt <- ex$features
  meta <- ex$meta
  ftype <- meta$ftype_levels[dt$ftype]
  lanes <- meta$feature_lanes
  lane_of <- setNames(rep(names(lanes), lengths(lanes)), unlist(lanes, use.names = FALSE))
  group <- unname(lane_of[ftype])
  lane_i <- match(group, names(lanes))
  nl <- length(lanes)
  aa <- ifelse(dt$s == dt$e, paste0("residue ", dt$s), paste0("residues ", dt$s, "-", dt$e))
  hover <- sprintf("%s<br>%s: <b>%s</b>%s<br>%s", dt$label, group, ftype,
                   ifelse(nzchar(dt$note), paste0(" (", dt$note, ")"), ""), aa)
  traces <- explain_bar_traces(dt, group, explain_feature_lane_colors(), rep(0.8, nrow(dt)),
                               nl - lane_i + 0.1, hover)
  explain_panel(traces, display_range, "UniProt", c(0, nl))
}

#' Hover text of translon calls.
#' @noRd
explain_translons_hover <- function(dt, meta) {
  type <- meta$type_levels[dt$type]
  tier <- meta$tier_levels[dt$tier]
  codon <- meta$codon_levels[dt$start_codon]
  len <- dt$end - dt$start + 1L
  num <- function(x, fmt) ifelse(is.na(x), "NA", sprintf(fmt, x))
  cons <- ifelse(is.na(dt$conserved), "NA", ifelse(dt$conserved, "yes", "no"))
  # meta$conserved_label: hover label of the conservation column (default: human lane, conserved in mouse);
  # NA / "" = the lane has no conservation column, line omitted
  cl <- if (is.null(meta$conserved_label)) "conserved translon (mouse)" else meta$conserved_label
  cons <- if (length(cl) != 1 || is.na(cl) || !nzchar(cl)) rep("", nrow(dt)) else paste0("<br>", cl, ": ", cons)
  eff <- if (!is.null(dt$rt_eff)) dt$rt_eff else rep(NA_real_, nrow(dt))
  genomic <- if (!is.null(dt$genomic)) dt$genomic else rep(NA_character_, nrow(dt))
  sprintf("%s<br><b>%s</b> (%s tier)<br>q = %s, P = %s<br>start codon %s; tx nt %d-%d (%d nt)%s%s%s",
          dt$label, type, tier, num(dt$q, "%.3g"), num(dt$P, "%.3g"), ifelse(is.na(codon), "NA", codon),
          dt$start, dt$end, len, cons,
          ifelse(is.na(eff), "", sprintf("<br>readthrough efficiency: %.3g", eff)),
          ifelse(is.na(genomic), "", paste0("<br>", genomic)))
}

#' Translon calls panel: one block of sub-lanes per transcript, boxes coloured
#' by type, opacity by tier.
#' @noRd
explain_translons_panel <- function(ex, display_range) {
  dt <- ex$translons
  meta <- ex$translons_meta
  txs <- unique(dt$tx)
  nsl <- dt$n_slots[match(txs, dt$tx)]
  offset <- setNames(c(0L, cumsum(nsl)[-length(nsl)]), txs)
  total <- sum(nsl)
  row <- offset[dt$tx] + dt$slot          # 1 = top row
  base <- total - row + 0.1
  type <- meta$type_levels[dt$type]
  tier <- meta$tier_levels[dt$tier]
  op <- meta$tier_opacity
  if (is.null(op)) op <- c(high = 0.95, candidate = 0.4)
  opacity <- unname(op[tier])
  opacity[is.na(opacity)] <- 0.4
  hover <- explain_translons_hover(dt, meta)
  cols <- meta$type_colors
  miss <- setdiff(unique(type), names(cols))
  if (length(miss)) cols[miss] <- "#B15928"
  traces <- explain_bar_traces(dt, type, cols, rep(0.8, nrow(dt)), base, hover)
  for (k in seq_along(traces)) {
    sel <- which(type == traces[[k]]$name)
    traces[[k]]$marker$opacity <- I(opacity[sel])
    traces[[k]]$marker$line <- list(width = 0.5, color = unname(cols[traces[[k]]$name]))
    traces[[k]]$meta <- list(rc_explain = "translons")
  }
  explain_panel(traces, display_range, "ORFs", c(0, total))  # short: the label sits in the ~60 px left margin
}

#' Translon type colour legend (HTML), NULL without a translon lane.
#' @noRd
explain_translons_legend_ui <- function(meta) {
  if (is.null(meta)) return(NULL)
  cols <- meta$type_colors
  items <- lapply(names(cols), function(s) tags$span(
    style = "display:inline-block; margin-right:10px; white-space:nowrap;",
    tags$span(style = paste0("display:inline-block; width:12px; height:12px; margin-right:4px; ",
                             "vertical-align:middle; background:", cols[[s]], ";")), s))
  tags$div(class = "explain-translons-legend", style = "font-size:12px;",
           tags$b(paste0("Translon calls (", meta$version, "; solid = high, faint = candidate): ")), items)
}

#' Bottom panels of the explainability mode (named explain_*).
#' @noRd
explain_bottom_panels <- function(ex, display_range) {
  if (is.null(ex)) return(list())
  panels <- list()
  if (nrow(ex$status)) panels$explain_status <- explain_status_panel(ex, display_range)
  if (nrow(ex$plddt)) panels$explain_plddt <- explain_plddt_panel(ex, display_range)
  if (nrow(ex$features)) panels$explain_features <- explain_features_panel(ex, display_range)
  if (!is.null(ex$translons) && nrow(ex$translons)) {
    panels$explain_translons <- explain_translons_panel(ex, display_range)
  }
  if (!is.null(ex$paraloc) && nrow(ex$paraloc)) {
    panels$explain_paraloc <- paraloc_panel(ex, display_range)
  }
  panels
}

#' Axis settings for explain panels (keeps hover, unlike other custom tracks).
#' @noRd
automateTicksExplainTrack <- function(plot) {
  plot  # explain_panel() already sets the axes; plotly::layout() would force a rebuild
}

# ---- expected coverage overlay --------------------------------------------

#' Scaled expected-coverage step line for one library profile.
#'
#' E is scaled so that its total over the displayed windows equals the library's
#' coverage at the usable positions of those windows.
#' @return list(x, y) or NULL
#' @noRd
explain_expected_line <- function(expected, profile, log_scale = FALSE) {
  if (is.null(expected) || !nrow(expected) || is.null(profile) || !nrow(profile)) return(NULL)
  count <- profile$count
  if (log_scale) count <- 2^count - 1
  pos <- as.integer(profile$position)
  wmax <- max(c(pos, expected$x1))
  v <- numeric(wmax)
  ok <- !is.na(pos) & pos >= 1L
  agg <- rowsum(count[ok], pos[ok])
  v[as.integer(rownames(agg))] <- agg[, 1]
  cv <- c(0, cumsum(v))
  observed <- sum(cv[expected$x1 + 1L] - cv[expected$x0])
  expected_total <- expected$E_shown_total[1]
  if (!is.finite(observed) || observed <= 0 || !is.finite(expected_total) || expected_total <= 0) return(NULL)
  y <- expected$density * observed / expected_total
  if (log_scale) y <- log2(y + 1)
  n <- nrow(expected)
  list(x = as.vector(rbind(expected$x0 - 0.5, expected$x1 + 0.5, NA)),
       y = as.vector(rbind(y, y, NA)))
}

#' Append the expected line to each Ribo-seq library track plot.
#' @noRd
explain_add_expected_traces <- function(plots, expected, profiles, withFrames,
                                        log_scale = FALSE, summary_track = FALSE,
                                        frames_type = "lines") {
  if (is.null(expected) || !nrow(expected) || identical(frames_type, "animate")) return(plots)
  nlib <- length(profiles)
  for (k in seq_len(nlib)) {
    if (!isTRUE(withFrames[k])) next
    i <- if (isTRUE(summary_track)) length(plots) - k + 1L else k
    line <- explain_expected_line(expected, profiles[[k]], isTRUE(rep(log_scale, length.out = nlib)[k]))
    if (is.null(line) || i < 1L || i > length(plots)) next
    p <- plots[[i]]
    if (is.null(p$x$data) || length(p$x$layoutAttrs)) p <- suppressWarnings(plotly::plotly_build(p))
    p$x$data[[length(p$x$data) + 1L]] <- list(
      type = "scatter", mode = "lines", x = I(line$x), y = I(line$y), xaxis = "x", yaxis = "y",
      line = list(color = "#7b3294", width = 2, shape = "hv"),
      connectgaps = FALSE, name = "Expected (E)", hoverinfo = "skip", showlegend = FALSE,
      meta = list(rc_explain = "expected"))
    plots[[i]] <- p
  }
  plots
}

# ---- controller -----------------------------------------------------------

#' Controller snapshot for explainability mode (plain list, cache safe).
#' @noRd
explain_controller <- function(input, dff, tx_grl, display_region) {
  off <- list(active = FALSE)
  if (!isTRUE(browser_input_or_default(input, "explain_mode", FALSE))) return(off)
  if (is.null(explain_store_dir())) return(off)
  org <- tryCatch(organism(dff), error = function(e) NA_character_)
  store <- tryCatch(explain_store_for(org), error = function(e) NULL)
  if (is.null(store)) {
    if (!is.null(shiny::getDefaultReactiveDomain())) {
      shiny::showNotification(paste0("Explainability mode is only available for ",
                                     explain_organism_labels(), "."), type = "warning")
    }
    return(off)
  }
  layers <- browser_input_or_default(input, "explain_layers", c("status", "expected", "exons"))
  viewMode <- isTRUE(browser_input_or_default(input, "viewMode", FALSE))
  other_tx <- isTRUE(browser_input_or_default(input, "other_tx", FALSE))
  ids <- explain_candidate_ids(store, shiny::isolate(input$tx), tx_grl, display_region,
                               include_overlapping = viewMode && other_tx, layers = layers)
  if (!length(ids)) {
    if (!is.null(shiny::getDefaultReactiveDomain())) {
      shiny::showNotification("Explainability: no transcript with explainer results at this locus.",
                              type = "message")
    }
  }
  list(active = length(ids) > 0, layers = layers, ids = ids, store_dir = store$dir,
       tx_exons = if (length(ids)) tx_grl[ids] else NULL)
}

#' Hash part for explainability inputs.
#' @noRd
explain_hash <- function(input) {
  paste(isTRUE(browser_input_or_default(input, "explain_mode", FALSE)),
        paste(browser_input_or_default(input, "explain_layers", ""), collapse = ","),
        sep = ":")
}

#' Settings UI (empty when no store is configured).
#' @noRd
explain_settings_ui <- function(ns) {
  if (is.null(explain_store_dir())) return(NULL)
  stores <- explain_stores_open()
  tl <- explain_stores_translons(stores)
  fluidRow(
    column(4, checkboxInput(ns("explain_mode"),
                            paste0("Explainability (", explain_organism_labels(names(stores)), ")"), FALSE)),
    column(8, conditionalPanel(
      sprintf("input['%s'] == true", ns("explain_mode")),
      checkboxGroupInput(ns("explain_layers"), NULL, choices = explain_layer_choices(tl),
                         selected = c("status", "expected", "exons"), inline = TRUE),
      explain_status_legend_store(stores),
      explain_translons_legend_ui(tl$translons$meta),
      if (isTRUE(tl$has_paraloc)) paraloc_legend_ui(Filter(Negate(is.null), lapply(stores, function(s) s$paraloc))[[1]]$meta)))
  )
}

#' Pseudo store for the settings UI: translon lane meta merged over all stores
#' (type colours union, version "v4_3 human, v1 mouse"); list(translons = NULL) without lanes.
#' @noRd
explain_stores_translons <- function(stores) {
  has_paraloc <- any(vapply(stores, function(s) !is.null(s$paraloc), logical(1)))
  metas <- lapply(stores, function(s) s$translons$meta)
  metas <- metas[!vapply(metas, is.null, logical(1))]
  if (!length(metas)) return(list(translons = NULL, has_paraloc = has_paraloc))
  if (length(metas) == 1) return(list(translons = list(meta = metas[[1]]), has_paraloc = has_paraloc))
  cols <- character()
  for (m in metas) cols <- c(cols, m$type_colors[setdiff(names(m$type_colors), names(cols))])
  ver <- paste(vapply(names(metas), function(o) paste(metas[[o]]$version, explain_organism_labels(o)), ""),
               collapse = ", ")
  list(translons = list(meta = list(version = ver, type_colors = cols)), has_paraloc = has_paraloc)
}

#' Status legend for the configured stores (union of their status levels;
#' NULL when none can be read).
#' @noRd
explain_status_legend_store <- function(stores = explain_stores_open()) {
  if (!length(stores)) return(NULL)
  meta <- stores[[1]]$meta
  for (s in stores[-1]) {
    m <- s$meta
    meta$status_levels <- union(meta$status_levels, m$status_levels)
    meta$status_colors <- c(meta$status_colors, m$status_colors[setdiff(names(m$status_colors), names(meta$status_colors))])
    meta$status_labels <- c(meta$status_labels, m$status_labels[setdiff(names(m$status_labels), names(meta$status_labels))])
    meta$type_colors <- c(meta$type_colors, m$type_colors[setdiff(names(m$type_colors), names(meta$type_colors))])
    meta$type_labels <- c(meta$type_labels, m$type_labels[setdiff(names(m$type_labels), names(meta$type_labels))])
  }
  explain_status_legend_ui(meta)
}

#' Enable the explainability switch only for supported organisms.
#' @noRd
explain_organism_observer <- function(input, session, df) {
  if (is.null(explain_store_dir())) return(invisible(NULL))
  observe({
    org <- tryCatch(organism(df()), error = function(e) NA_character_)
    ok <- explain_organism_supported(org)
    shinyjs::toggleState("explain_mode", condition = ok)
    if (!ok && isTRUE(isolate(input$explain_mode))) {
      updateCheckboxInput(session, "explain_mode", value = FALSE)
    }
  })
}
