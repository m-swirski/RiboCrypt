# Explainability mode for the single-locus browser.
#
# Shows precomputed explainer layers for MANE transcripts (human for now):
#   - window status ribbon (90-nt windows, status_fam, hover: type, O, E, log2 O/E)
#   - expected P-site coverage E overlaid on Ribo-seq tracks (scaled to each library)
#   - AlphaFold pLDDT per codon
#   - UniProt features per protein (domains, disorder, topology, sites, AF helix/strand)
# Data come from a compact store directory (built at release time by
# build_explain_store.R): index.fst (one row per transcript with row ranges),
# windows.fst, lowmap.fst, plddt.fst, features.fst and meta.rds. Rows are read
# lazily per transcript and kept in an in-memory cache.
#
# The store directory is given by option "RiboCrypt.explain_store" or the
# environment variable RIBOCRYPT_EXPLAIN_STORE; without it the mode is hidden.

#' @noRd
explain_layer_choices <- function() {
  c("Window status" = "status", "Expected E on P-sites" = "expected",
    "AlphaFold pLDDT" = "plddt", "UniProt features" = "features")
}

#' Explainability store directory, NULL when not configured.
#' @noRd
explain_store_dir <- function() {
  d <- getOption("RiboCrypt.explain_store", Sys.getenv("RIBOCRYPT_EXPLAIN_STORE", ""))
  if (is.null(d) || !nzchar(d) || !dir.exists(d) ||
      !file.exists(file.path(d, "index.fst"))) return(NULL)
  normalizePath(d)
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
                cache = new.env(parent = emptyenv()))
  assign(key, store, envir = .rc_explain_env)
  store
}

#' Is the organism supported by the store?
#' @noRd
explain_organism_supported <- function(organism_name, store = explain_store_open()) {
  if (is.null(store) || length(organism_name) != 1 || is.na(organism_name)) return(FALSE)
  norm <- function(x) tolower(gsub("[_ ]+", " ", trimws(x)))
  identical(norm(organism_name), norm(store$meta$organism))
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
                                  include_overlapping = FALSE, max_tx = 5L) {
  if (is.null(store)) return(character())
  ids <- intersect(selected, names(store$row))
  if (include_overlapping) {
    spans <- explain_store_spans(store, tx_grl)
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
explain_store_spans <- function(store, tx_grl) {
  key <- paste0("spans:", length(tx_grl), ":", names(tx_grl)[1], ":",
                names(tx_grl)[length(tx_grl)])
  spans <- store$cache[[key]]
  if (is.null(spans)) {
    ids <- intersect(names(store$row), names(tx_grl))
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
explain_locus_data <- function(explain, display_range, store = explain_store_open()) {
  if (is.null(store) || !isTRUE(explain$active) || !length(explain$ids)) return(NULL)
  layers <- explain$layers
  display <- unlist(display_range, use.names = FALSE)
  out <- list(status = list(), expected = list(), plddt = list(), features = list(),
              meta = store$meta, ids = character())
  for (lane in seq_along(explain$ids)) {
    id <- explain$ids[lane]
    d <- explain_store_tx(store, id)
    if (is.null(d) || !(id %in% names(explain$tx_exons))) next
    exons <- explain$tx_exons[[id]]
    if (sum(width(exons)) != d$tx$n) {
      warning("Explainability: transcript length mismatch for ", id, ", skipped")
      next
    }
    map <- explain_tx_display_map(exons, display)
    if (!nrow(map)) next
    out$ids <- c(out$ids, id)
    label <- paste0(d$tx$gene, " (", id, ")")
    if (any(c("status", "expected") %in% layers) && !is.null(d$windows)) {
      w <- d$windows
      if ("status" %in% layers) {
        p <- explain_map_ranges(w$start, w$end, map)
        if (nrow(p)) out$status[[id]] <- data.table(p, w[p$id], lane = lane, label = label)
      }
      if ("expected" %in% layers) out$expected[[id]] <- explain_expected_segments(d, map)
    }
    if ("plddt" %in% layers && !is.null(d$plddt)) {
      out$plddt[[id]] <- explain_plddt_bins(d, map, lane = lane, label = label)
    }
    if ("features" %in% layers && !is.null(d$features)) {
      f <- d$features
      cs <- d$tx$cs
      p <- explain_map_ranges(cs + 3L * (f$s - 1L), cs + 3L * f$e - 1L, map)
      if (nrow(p)) out$features[[id]] <- data.table(p, f[p$id], label = label)
    }
  }
  if (!length(out$ids)) return(NULL)
  out$status <- rbindlist(out$status)
  out$expected <- rbindlist(out$expected)
  out$plddt <- rbindlist(out$plddt)
  out$features <- rbindlist(out$features)
  out
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
  yaxis <- custom_seq_track_yaxis(ylab)
  yaxis$range <- yrange
  yaxis$autorange <- FALSE
  yaxis$title$font <- list(size = 22)
  xaxis <- custom_seq_track_xaxis(display_range)
  xaxis$autorange <- FALSE
  p$x$data <- traces
  p$x$layout$xaxis <- xaxis
  p$x$layout$yaxis <- yaxis
  p$x$layout$bargap <- 0
  p$x$layout$barmode <- "overlay"
  p$x$layout$margin <- list(t = 0, r = 0, b = 0, l = 0, pad = 0)
  p$x$layout$showlegend <- FALSE
  p
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
  hover <- sprintf(
    "%s<br>%s window %d-%d nt: <b>%s</b>%s<br>O = %s, E = %s, log2 O/E = %.2f<br>all-mapper track: %s",
    dt$label, region, dt$start, dt$end, status,
    ifelse(is.na(type), "", paste0("<br>", type)),
    format(round(dt$O), big.mark = ",", trim = TRUE), format(round(dt$E), big.mark = ",", trim = TRUE),
    log2((dt$O + 1) / (dt$E + 1)), ifelse(is.na(status_all), "NA", status_all))
  nl <- max(dt$lane)
  traces <- explain_bar_traces(dt, status, meta$status_colors, rep(0.9, nrow(dt)),
                               nl - dt$lane + 0.05, hover)
  explain_panel(traces, display_range, "Explain", c(0, nl))
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

#' Bottom panels of the explainability mode (named explain_*).
#' @noRd
explain_bottom_panels <- function(ex, display_range) {
  if (is.null(ex)) return(list())
  panels <- list()
  if (nrow(ex$status)) panels$explain_status <- explain_status_panel(ex, display_range)
  if (nrow(ex$plddt)) panels$explain_plddt <- explain_plddt_panel(ex, display_range)
  if (nrow(ex$features)) panels$explain_features <- explain_features_panel(ex, display_range)
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
  store <- explain_store_open()
  if (is.null(store)) return(off)
  org <- tryCatch(organism(dff), error = function(e) NA_character_)
  if (!explain_organism_supported(org, store)) {
    if (!is.null(shiny::getDefaultReactiveDomain())) {
      shiny::showNotification(paste0("Explainability mode is only available for ",
                                     store$meta$organism, " (MANE transcripts)."), type = "warning")
    }
    return(off)
  }
  layers <- browser_input_or_default(input, "explain_layers", c("status", "expected"))
  viewMode <- isTRUE(browser_input_or_default(input, "viewMode", FALSE))
  other_tx <- isTRUE(browser_input_or_default(input, "other_tx", FALSE))
  ids <- explain_candidate_ids(store, shiny::isolate(input$tx), tx_grl, display_region,
                               include_overlapping = viewMode && other_tx)
  if (!length(ids)) {
    if (!is.null(shiny::getDefaultReactiveDomain())) {
      shiny::showNotification("Explainability: no MANE transcript with explainer results at this locus.",
                              type = "message")
    }
  }
  list(active = length(ids) > 0, layers = layers, ids = ids,
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
  fluidRow(
    column(4, checkboxInput(ns("explain_mode"), "Explainability (MANE, human)", FALSE)),
    column(8, conditionalPanel(
      sprintf("input['%s'] == true", ns("explain_mode")),
      checkboxGroupInput(ns("explain_layers"), NULL, choices = explain_layer_choices(),
                         selected = c("status", "expected"), inline = TRUE)))
  )
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
