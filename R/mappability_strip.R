# Mappability strip: browser bottom track drawing 28-mer mappability as coloured blocks.
#
# Two rows from a compact store (rc_mappability_v*, built by build_mappability_store_v1.R):
#   "unique (perfect match)": exact 28-mer uniqueness (genome + exon-junction windows; canonical transcripts)
#   "<=2 mismatches":         genome-wide 28-mer track, unique = no other site within 2 mismatches
# Values are P-site positions of a 28-nt read (offset 12). Colour = fraction of uniquely mappable
# P-site positions per block (sequential palette, light = poorly mappable); positions outside the
# store domain (introns, intergenic, non-canonical transcripts for the exact row) are left empty.
#
# Store lookup: option "RiboCrypt.mappability_store" (directory, or named vector by organism), else
# the newest <reference dir>/mappability/rc_mappability_v* next to the experiment's genome fasta.
# Without a store the browser keeps the old bigWig bar track (28mers_mappability.bw).

.rc_mappability_env <- new.env(parent = emptyenv())

#' Mappability store directory for an experiment, NULL when none exists.
#' @noRd
mappability_store_dir <- function(df) {
  opt <- getOption("RiboCrypt.mappability_store", NULL)
  if (!is.null(opt)) {
    if (!is.null(names(opt)) && any(nzchar(names(opt)))) {
      org <- tryCatch(organism(df), error = function(e) NA_character_)
      hit <- opt[explain_norm_organism(names(opt)) == explain_norm_organism(org)]
      opt <- if (length(hit)) hit[[1]] else NULL
    } else opt <- opt[[1]]
    if (!is.null(opt) && file.exists(file.path(opt, "meta.rds"))) return(normalizePath(opt))
    return(NULL)
  }
  fa <- tryCatch(df@fafile, error = function(e) NULL)
  if (is.null(fa) || !length(fa) || !nzchar(fa)) return(NULL)
  mappability_find_store(file.path(dirname(fa), "mappability"))
}

#' Newest rc_mappability_v<N> directory with a meta.rds in a folder, NULL when none.
#' @noRd
mappability_find_store <- function(dir) {
  if (!dir.exists(dir)) return(NULL)
  d <- list.files(dir, pattern = "^rc_mappability_v[0-9]+$", full.names = TRUE)
  d <- d[file.exists(file.path(d, "meta.rds")) & file.exists(file.path(d, "index.fst"))]
  if (!length(d)) return(NULL)
  v <- as.integer(sub("^rc_mappability_v", "", basename(d)))
  normalizePath(d[which.max(v)])
}

#' Open (and cache) a mappability store: meta + bin index.
#' @noRd
mappability_store_open <- function(dir) {
  if (is.null(dir)) return(NULL)
  key <- paste0("store:", dir)
  s <- .rc_mappability_env[[key]]
  if (!is.null(s)) return(s)
  s <- tryCatch(list(dir = dir, meta = readRDS(file.path(dir, "meta.rds")),
                     index = fst::read_fst(file.path(dir, "index.fst"), as.data.table = TRUE)),
                error = function(e) { warning("Mappability store unreadable: ", conditionMessage(e)); NULL })
  if (!is.null(s)) assign(key, s, envir = .rc_mappability_env)
  s
}

#' Match display seqnames to store seqlevels ("chr1" / "1", "chrM" / "MT").
#' @noRd
mappability_chr_code <- function(chr, seqlevels) {
  i <- match(chr, seqlevels)
  if (is.na(i)) i <- match(sub("^chr", "", chr), seqlevels)
  if (is.na(i) && chr %in% c("chrM", "M")) i <- match("MT", seqlevels)
  if (is.na(i) && chr == "MT") i <- match("chrM", seqlevels)
  i
}

#' Rows of one store table (track x part) overlapping a genomic range of one chromosome.
#' @return data.table(start, end, strand[, n]) or NULL
#' @noRd
mappability_read <- function(store, track, part, chr_code, from, to) {
  bin <- store$meta$bin
  b <- seq.int((max(1L, from) - 1L) %/% bin, (to - 1L) %/% bin)
  ix <- store$index[store$index$track == track & store$index$part == part &
                      store$index$chr == chr_code & store$index$bin %in% b]
  if (!nrow(ix)) return(NULL)
  ix <- ix[order(ix$from)]
  # contiguous bins -> one read each
  brk <- c(TRUE, ix$from[-1L] != ix$to[-nrow(ix)] + 1L)
  g <- cumsum(brk)
  file <- file.path(store$dir, paste0(track, "_", part, ".fst"))
  res <- lapply(split(seq_len(nrow(ix)), g), function(k) {
    fst::read_fst(file, from = ix$from[k[1]], to = ix$to[k[length(k)]], as.data.table = TRUE)
  })
  d <- data.table::rbindlist(res)
  d[d$end >= from & d$start <= to]
}

#' Map genomic intervals to browser display x (1-based position along the display region).
#' @param display GRanges, display region in display order (one strand)
#' @return data.table(id, x0, x1) with x0 <= x1, split at display segment borders
#' @noRd
mappability_genomic_to_x <- function(start, end, display) {
  empty <- data.table(id = integer(), x0 = integer(), x1 = integer())
  if (!length(start) || !length(display)) return(empty)
  ds <- start(display); de <- end(display)
  dcum <- cumsum(c(0L, (de - ds + 1L)[-length(ds)]))
  od <- order(ds)
  oq <- order(start)
  hit <- explain_intersect_sorted(as.integer(start[oq]), as.integer(end[oq]), ds[od], de[od])
  if (!length(hit$qi)) return(empty)
  q <- oq[hit$qi]; s <- od[hit$si]
  minus <- as.character(strand(display))[1] == "-"
  if (minus) {
    x0 <- dcum[s] + de[s] - hit$hi + 1L; x1 <- dcum[s] + de[s] - hit$lo + 1L
  } else {
    x0 <- dcum[s] + hit$lo - ds[s] + 1L; x1 <- dcum[s] + hit$hi - ds[s] + 1L
  }
  data.table(id = q, x0 = as.integer(x0), x1 = as.integer(x1))
}

#' Per-x coverage (0/1) of mapped intervals on 1..L via a difference array.
#' @noRd
mappability_cover <- function(x, L) {
  v <- integer(L + 1L)
  if (nrow(x)) {
    a <- tabulate(x$x0, nbins = L + 1L); b <- tabulate(x$x1 + 1L, nbins = L + 1L)
    v <- cumsum(a - b)
  }
  as.integer(v[seq_len(L)] > 0L)
}

#' Block summaries of one track for the display region.
#'
#' @param max_blocks number of blocks at most (one block = ceiling(L / max_blocks) nt)
#' @return data.table(x0, x1, frac, defined, n_max) for blocks with defined positions, or NULL
#' @noRd
mappability_track_blocks <- function(store, track, display, max_blocks = 1500L) {
  tr <- store$meta$tracks[[track]]
  if (is.null(tr)) return(NULL)
  L <- sum(width(display))
  if (L < 1L) return(NULL)
  minus <- as.character(strand(display))[1] == "-"
  shift <- if (minus) tr$shift[["minus"]] else tr$shift[["plus"]]
  st_code <- if (minus) 2L else 1L
  chr <- as.character(seqnames(display))[1]
  cc <- mappability_chr_code(chr, store$meta$seqlevels)
  if (is.na(cc)) return(NULL)
  from <- min(start(display)) - shift; to <- max(end(display)) - shift
  get <- function(part) {
    d <- mappability_read(store, track, part, cc, from, to)
    if (is.null(d) || !nrow(d)) return(NULL)
    d <- d[d$strand == 0L | d$strand == st_code]
    if (shift != 0L) { d$start <- d$start + shift; d$end <- d$end + shift }
    d
  }
  def <- get("def")
  if (is.null(def)) return(NULL)
  nu <- get("nu")
  xd <- mappability_genomic_to_x(def$start, def$end, display)
  if (!nrow(xd)) return(NULL)
  cov_def <- mappability_cover(xd, L)
  cov_nu <- integer(L); nmax <- integer(L)
  if (!is.null(nu) && nrow(nu)) {
    xn <- mappability_genomic_to_x(nu$start, nu$end, display)
    cov_nu <- mappability_cover(xn, L) * cov_def
    if (!is.null(nu$n) && nrow(xn) && any(!is.na(nu$n))) {
      # max sites per position: only needed per block, so take it per interval
      nn <- nu$n[xn$id]; nn[is.na(nn)] <- 2L
      nmax_iv <- data.table(x0 = xn$x0, x1 = xn$x1, n = nn)
    } else nmax_iv <- NULL
  } else nmax_iv <- NULL
  bs <- max(1L, as.integer(ceiling(L / max_blocks)))
  blk <- (seq_len(L) - 1L) %/% bs
  nb <- blk[L] + 1L
  def_b <- tabulate(blk[cov_def == 1L] + 1L, nbins = nb)
  uniq_b <- tabulate(blk[cov_def == 1L & cov_nu == 0L] + 1L, nbins = nb)
  b <- which(def_b > 0L)
  if (!length(b)) return(NULL)
  x0 <- (b - 1L) * bs + 1L; x1 <- pmin(b * bs, L)
  res <- data.table(x0 = x0, x1 = x1, frac = uniq_b[b] / def_b[b], defined = def_b[b] / (x1 - x0 + 1L),
                    n_max = NA_integer_)
  if (!is.null(nmax_iv)) {
    # largest N of any non-unique interval touching the block
    b0 <- (nmax_iv$x0 - 1L) %/% bs + 1L; b1 <- (nmax_iv$x1 - 1L) %/% bs + 1L
    k <- rep.int(seq_len(nrow(nmax_iv)), b1 - b0 + 1L)
    bb <- sequence(b1 - b0 + 1L, from = b0)
    mx <- tapply(nmax_iv$n[k], bb, max)
    res$n_max <- as.integer(mx[match(b, as.integer(names(mx)))])
  }
  res
}

#' Fraction-unique classes and colours (sequential, light = poorly mappable).
#' @noRd
mappability_classes <- function() {
  c("< 20% unique" = "#ffffcc", "20-40%" = "#c7e9b4", "40-60%" = "#7fcdbb",
    "60-80%" = "#41b6c4", "80-99%" = "#2c7fb8", "100% unique" = "#253494")
}
mappability_class_of <- function(frac) {
  lab <- names(mappability_classes())
  lab[findInterval(frac, c(0.2, 0.4, 0.6, 0.8, 1 - 1e-9)) + 1L]
}

#' Merge adjacent blocks of the same class into one segment.
#' @noRd
mappability_merge_blocks <- function(b) {
  if (is.null(b) || !nrow(b)) return(b)
  cls <- mappability_class_of(b$frac)
  n <- nrow(b)
  new <- c(TRUE, cls[-1L] != cls[-n] | b$x0[-1L] != b$x1[-n] + 1L)
  g <- cumsum(new)
  w <- (b$x1 - b$x0 + 1) * b$defined
  data.table(x0 = as.integer(tapply(b$x0, g, min)), x1 = as.integer(tapply(b$x1, g, max)),
             frac = as.numeric(tapply(b$frac * w, g, sum) / tapply(w, g, sum)),
             defined = as.numeric(tapply(w, g, sum) / tapply(b$x1 - b$x0 + 1, g, sum)),
             n_max = suppressWarnings(as.integer(tapply(b$n_max, g, function(v) if (all(is.na(v))) NA else max(v, na.rm = TRUE)))),
             cls = cls[new])
}

#' Strip data for the display region: one segment table per store track.
#' @param display_range GRangesList (one group) or GRanges, display region in display order
#' @return list(tracks = named list of data.tables, meta) or NULL
#' @noRd
mappability_strip_data <- function(store, display_range, max_blocks = 1500L) {
  if (is.null(store)) return(NULL)
  display <- if (is(display_range, "GRangesList")) unlist(display_range, use.names = FALSE) else display_range
  if (!length(display)) return(NULL)
  tracks <- list()
  for (tr in intersect(c("exact", "le2mm"), names(store$meta$tracks))) {
    b <- tryCatch(mappability_track_blocks(store, tr, display, max_blocks), error = function(e) {
      warning("Mappability track ", tr, ": ", conditionMessage(e)); NULL })
    m <- mappability_merge_blocks(b)
    if (!is.null(m) && nrow(m)) tracks[[tr]] <- m
  }
  if (!length(tracks)) return(NULL)
  list(tracks = tracks, meta = store$meta)
}

#' Plotly panel of the mappability strip (two rows, coloured blocks).
#' @noRd
mappability_strip_panel <- function(sd, display_range) {
  rows <- intersect(c("exact", "le2mm"), names(sd$meta$tracks))   # top to bottom
  nr <- length(rows)
  cols <- mappability_classes()
  traces <- list()
  ann <- list()
  for (i in seq_along(rows)) {
    tr <- rows[i]
    info <- sd$meta$tracks[[tr]]
    base <- nr - i + 0.08
    ann[[length(ann) + 1L]] <- list(text = info$short, xref = "paper", x = 0, xanchor = "right", xshift = -4,
                                    yref = "y", y = base + 0.42, yanchor = "middle", showarrow = FALSE,
                                    font = list(size = 9))
    dt <- sd$tracks[[tr]]
    if (is.null(dt)) next
    ncl <- info$n_classes
    n_lab <- ifelse(dt$n_max >= 65535L, "\u2265 1000 (saturated)", as.character(dt$n_max))
    if (length(ncl)) {
      hit <- as.character(dt$n_max) %in% names(ncl)
      n_lab[hit] <- ncl[as.character(dt$n_max[hit])]
    }
    hover <- sprintf("%s 28-mer mappability<br>%.0f%% of P-site positions unique%s%s<br>display nt %d-%d",
                     info$label, 100 * dt$frac,
                     ifelse(dt$defined < 0.999, sprintf(" (%.0f%% of positions in store domain)", 100 * dt$defined), ""),
                     ifelse(is.na(dt$n_max) | tr == "exact", "", paste0("<br>max sites within 2 mm: ", n_lab)),
                     dt$x0, dt$x1)
    tt <- explain_bar_traces(dt, dt$cls, cols, rep(0.84, nrow(dt)), rep(base, nrow(dt)), hover)
    for (k in seq_along(tt)) tt[[k]]$meta <- list(rc_track = "mappability", rc_row = tr)
    traces <- c(traces, tt)
  }
  p <- explain_panel(traces, display_range, "", c(0, nr))
  p$x$layout$annotations <- ann
  p
}

#' Mappability bottom panel for an experiment, NULL without a store / data at the locus.
#' @noRd
mappability_strip_track_panel <- function(df, display_range) {
  store <- mappability_store_open(mappability_store_dir(df))
  if (is.null(store)) return(NULL)
  print("- Loading mappability strip")
  sd <- mappability_strip_data(store, display_range)
  if (is.null(sd)) return(NULL)
  mappability_strip_panel(sd, display_range)
}

#' Colour legend of the mappability strip (HTML).
#' @noRd
mappability_legend_ui <- function() {
  cols <- mappability_classes()
  items <- lapply(names(cols), function(s) tags$span(
    style = "display:inline-block; margin-right:8px; white-space:nowrap;",
    tags$span(style = paste0("display:inline-block; width:12px; height:12px; margin-right:3px; ",
                             "vertical-align:middle; border:1px solid #bbb; background:", cols[[s]], ";")), s))
  tags$div(class = "mappability-legend", style = "font-size:12px;",
           tags$b("Mappability (28-mer P-sites; rows: unique perfect match / \u22642 mismatches): "), items)
}
