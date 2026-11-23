# Paralloc localisation lane of the explainability mode.
#
# Segments over the signal-bearing span of each product (store paraloc_track_v*, built by
# build_paraloc_store_v1.R from Michal's Paralloc predictor outputs): N-terminal signal peptide (SP)
# or mitochondrial targeting sequence (MTS) of annotated CDSs and translons, NTE gain / masking of a
# targeting signal, and (faint) the predicted compartment of translons / NTEs without a signal.
# Hover: predictor family, per-tool scores, span policy, product.
#
# Store files: index.fst (transcript, n, cs, ce, s_from, s_to, gene), segments.fst (start, end, kind,
# cls, score, delta, span, ttype, id, tools, genomic; transcript nt, 1-based), meta.rds (levels,
# colours, kind opacity, policy). Lookup: option "RiboCrypt.paraloc_store" (directory, or named
# vector by organism), else the newest
# <reference dir>/predicted_translons/localization/paraloc_track_v* of the organism's explain store
# (<reference dir>/explainability/<version>).

#' Paralloc store directory for an explain store (NULL when none).
#' @noRd
paraloc_store_dir <- function(explain_dir, organism = NULL) {
  opt <- getOption("RiboCrypt.paraloc_store", NULL)
  if (!is.null(opt)) {
    if (!is.null(names(opt)) && any(nzchar(names(opt)))) {
      hit <- opt[explain_norm_organism(names(opt)) == explain_norm_organism(organism)]
      opt <- if (length(hit)) hit[[1]] else NULL
    } else opt <- opt[[1]]
    if (!is.null(opt) && paraloc_store_ok(opt)) return(normalizePath(opt))
    return(NULL)
  }
  if (is.null(explain_dir)) return(NULL)
  loc <- file.path(dirname(dirname(explain_dir)), "predicted_translons", "localization")
  if (!dir.exists(loc)) return(NULL)
  d <- list.files(loc, pattern = "^paraloc_track_v[0-9]+$", full.names = TRUE)
  d <- d[vapply(d, paraloc_store_ok, logical(1))]
  if (!length(d)) return(NULL)
  normalizePath(d[which.max(as.integer(sub("^paraloc_track_v", "", basename(d))))])
}

#' @noRd
paraloc_store_ok <- function(dir) {
  all(file.exists(file.path(dir, c("index.fst", "segments.fst", "meta.rds"))))
}

#' Open a Paralloc store (index + meta; segments read per transcript), NULL when absent / unreadable.
#' @noRd
paraloc_open <- function(dir) {
  if (is.null(dir)) return(NULL)
  tryCatch({
    index <- fst::read_fst(file.path(dir, "index.fst"), as.data.table = TRUE)
    index <- index[index$s_to >= index$s_from]
    meta <- readRDS(file.path(dir, "meta.rds"))
    list(dir = dir, file = file.path(dir, "segments.fst"), index = index, meta = meta,
         row = setNames(seq_len(nrow(index)), index$transcript), cache = new.env(parent = emptyenv()))
  }, error = function(e) {
    warning("Explainability: could not read Paralloc store ", dir, ": ", conditionMessage(e))
    NULL
  })
}

#' Paralloc segments of one transcript with packed sub-lanes (lazy, cached); NULL when none.
#' @noRd
paraloc_tx <- function(pl, id, max_cache = 2000L) {
  if (is.null(pl)) return(NULL)
  hit <- pl$cache[[id]]
  if (!is.null(hit)) return(hit)
  i <- pl$row[id]
  if (is.na(i)) return(NULL)
  r <- pl$index[i]
  seg <- fst::read_fst(pl$file, from = r$s_from, to = r$s_to, as.data.table = TRUE)
  seg$slot <- paraloc_pack(seg$start, seg$end)
  res <- list(tx = r, seg = seg)
  if (length(ls(pl$cache)) >= max_cache) rm(list = ls(pl$cache), envir = pl$cache)
  assign(id, res, envir = pl$cache)
  res
}

#' Greedy interval packing: sub-lane per segment so that segments in one lane do not overlap.
#' @noRd
paraloc_pack <- function(start, end) {
  n <- length(start)
  slot <- integer(n)
  if (!n) return(slot)
  o <- order(start, -end)
  last <- integer()
  for (k in o) {
    free <- which(last < start[k])
    if (length(free)) { slot[k] <- free[1]; last[free[1]] <- end[k] }
    else { last <- c(last, end[k]); slot[k] <- length(last) }
  }
  slot
}

#' Paralloc segments of one transcript in display coordinates (split at exon junctions in genomic view).
#' @noRd
paraloc_locus <- function(pl, id, exons, display) {
  d <- paraloc_tx(pl, id)
  if (is.null(d) || !nrow(d$seg)) return(NULL)
  if (sum(width(exons)) != d$tx$n) {
    warning("Explainability: transcript length mismatch for ", id, " (Paralloc lane), skipped")
    return(NULL)
  }
  map <- explain_tx_display_map(exons, display)
  if (!nrow(map)) return(NULL)
  p <- explain_map_ranges(d$seg$start, d$seg$end, map)
  if (!nrow(p)) return(NULL)
  gene <- if (!is.null(d$tx$gene) && !is.na(d$tx$gene)) d$tx$gene else ""
  data.table(x0 = p$x0, x1 = p$x1, d$seg[p$id], n_slots = max(d$seg$slot), tx = id,
             label = paste0(gene, " (", id, ")"))
}

#' Plain-language labels of the Paralloc kinds, classes and span policies.
#' @noRd
paraloc_labels <- function() {
  list(kind = c(cds_signal = "N-terminal signal of the annotated CDS",
                orf_signal = "N-terminal signal of the translon",
                nte_gain = "N-terminal extension gains a targeting signal",
                nte_masks = "N-terminal extension masks the CDS signal (dual localisation candidate)",
                nte_context = "N-terminal extension, no targeting change (predicted compartment)",
                orf_whole = "translon without targeting signal (predicted compartment)"),
       cls = c(SP = "signal peptide (SP)", MTS = "mitochondrial targeting sequence (MTS)"),
       span = c(uniprot = "UniProt span", default_len = "default span (median UniProt length)",
                extension = "extension", "extension+default_len" = "extension + default span",
                whole_orf = "whole ORF"),
       family = c(SP = "SignalP6 + TargetP2 (mean)", MTS = "TargetP2 + MitoFates + MitoProt (mean)"))
}

#' Hover text of Paralloc segments.
#' @noRd
paraloc_hover <- function(dt, meta) {
  lab <- paraloc_labels()
  kind <- meta$kind_levels[dt$kind]
  cls <- meta$cls_levels[dt$cls]
  span <- meta$span_levels[dt$span]
  ttype <- meta$ttype_levels[dt$ttype]
  cls_lab <- ifelse(cls %in% names(lab$cls), lab$cls[cls], gsub("_", " ", cls))
  kind_lab <- ifelse(kind %in% names(lab$kind), lab$kind[kind], kind)
  span_lab <- ifelse(span %in% names(lab$span), lab$span[span], span)
  fam <- ifelse(cls %in% names(lab$family), paste0(lab$family[cls], ": "), "")
  tools <- if (!is.null(dt$tools)) dt$tools else rep(NA_character_, nrow(dt))
  pred <- ifelse(is.na(tools) | !nzchar(tools), "", paste0("<br>predictor: ", fam, tools))
  delta <- ifelse(is.na(dt$delta), "", sprintf(" (NTE+CDS minus CDS: %+.2f)", dt$delta))
  genomic <- if (!is.null(dt$genomic)) dt$genomic else rep(NA_character_, nrow(dt))
  sprintf("%s<br><b>%s</b><br>%s<br>score %.2f%s%s<br>%s %s; tx nt %d-%d (%s)%s",
          dt$label, cls_lab, kind_lab, dt$score, delta, pred, ttype, dt$id, dt$start, dt$end, span_lab,
          ifelse(is.na(genomic) | !nzchar(genomic), "", paste0("<br>", genomic)))
}

#' Paralloc lane: one block of sub-lanes per transcript, colour = class, opacity = kind.
#' @noRd
paraloc_panel <- function(ex, display_range) {
  dt <- ex$paraloc
  meta <- ex$paraloc_meta
  txs <- unique(dt$tx)
  nsl <- dt$n_slots[match(txs, dt$tx)]
  offset <- setNames(c(0L, cumsum(nsl)[-length(nsl)]), txs)
  total <- sum(nsl)
  row <- offset[dt$tx] + dt$slot
  base <- total - row + 0.1
  cls <- meta$cls_levels[dt$cls]
  kind <- meta$kind_levels[dt$kind]
  cols <- meta$cls_colors
  miss <- setdiff(unique(cls), names(cols))
  if (length(miss)) cols[miss] <- "#969696"
  op <- meta$kind_alpha
  opacity <- if (length(op)) unname(op[kind]) else rep(1, nrow(dt))
  opacity[is.na(opacity)] <- 0.6
  hover <- paraloc_hover(dt, meta)
  traces <- explain_bar_traces(dt, cls, cols, rep(0.8, nrow(dt)), base, hover)
  for (k in seq_along(traces)) {
    sel <- which(cls == traces[[k]]$name)
    traces[[k]]$marker$opacity <- I(opacity[sel])
    traces[[k]]$marker$line <- list(width = 0.5, color = unname(cols[traces[[k]]$name]))
    traces[[k]]$meta <- list(rc_explain = "paraloc")
  }
  explain_panel(traces, display_range, "Loc", c(0, total))
}

#' Paralloc class colour legend (HTML), NULL without meta.
#' @noRd
paraloc_legend_ui <- function(meta) {
  if (is.null(meta)) return(NULL)
  cols <- meta$cls_colors
  lab <- paraloc_labels()$cls
  items <- lapply(names(cols), function(s) tags$span(
    style = "display:inline-block; margin-right:10px; white-space:nowrap;",
    tags$span(style = paste0("display:inline-block; width:12px; height:12px; margin-right:4px; ",
                             "vertical-align:middle; background:", cols[[s]], ";")),
    if (s %in% names(lab)) lab[[s]] else gsub("_", " ", s)))
  tags$div(class = "explain-paraloc-legend", style = "font-size:12px;",
           tags$b(paste0("Localisation (Paralloc, ", meta$version, "; solid = targeting signal, faint = compartment only): ")),
           items)
}
