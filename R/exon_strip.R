# Exon strip: thin browser lane where the exons of the selected transcript alternate between two
# colours (odd exons colour 1, even exons colour 2), so exon junctions and empty exons are easy to
# see. Exons are numbered in transcript order (exon 1 = 5'-most exon of the transcript, on both
# strands). Transcript view: blocks over transcript coordinates, touching at the junctions.
# Genomic view: one block per exon, introns empty; with collapsed introns the gaps between exons
# are bridged by a thin connector so the strip reads as contiguous.
#
# Shown by the browser setting "Exons (alternating)" (input exon_strip) and by the explainability
# layer "exons" (placed directly below the window status lane).

#' Exon strip colours (odd / even exon).
#'
#' Two neutrals that do not compete with the frame colours (red / green / blue, colour-blind
#' theme) or the explain status palette: dark slate and warm stone. Validated with the dataviz
#' palette validator (light surface): lightness band PASS, adjacent CVD dE 22.3 (protan) PASS,
#' normal-vision dE 23.5 PASS; chroma floor fails by design (neutral). The stone is 2.7:1 against
#' white, so exon numbers are drawn on wide blocks as secondary encoding.
#' @noRd
exon_strip_colors <- function() c("odd exon" = "#4d5a6a", "even exon" = "#a29b90")

#' Text colour of exon labels on each block colour.
#' @noRd
exon_strip_label_colors <- function() c("odd exon" = "#ffffff", "even exon" = "#1f1f1f")

#' Exons of one transcript in transcript order (5' -> 3').
#' @noRd
exon_strip_tx_order <- function(exons) {
  if (!length(exons)) return(exons)
  minus <- as.character(strand(exons))[1] == "-"
  exons[order(start(exons), decreasing = minus)]
}

#' Exon (and collapsed-intron) segments of transcripts in display coordinates.
#'
#' @param tx_exons GRangesList of transcripts (names = transcript ids)
#' @param display_range GRangesList (one group) or GRanges, display region in display order
#' @param collapsed logical, introns collapsed (genomic view): add intron connector rows
#' @return list(exons = data.table(x0, x1, lane, tx, exon, n_exons, parity, seqnames, start, end,
#'   strand, width, tx_start, tx_end, shown), introns = data.table(x0, x1, lane, tx, intron, start,
#'   end, width)) or NULL when nothing is in view
#' @noRd
exon_strip_data <- function(tx_exons, display_range, collapsed = FALSE) {
  display <- if (is(display_range, "GRangesList")) unlist(display_range, use.names = FALSE) else display_range
  if (!length(tx_exons) || !length(display)) return(NULL)
  ex_out <- list(); in_out <- list()
  lane <- 0L
  for (id in names(tx_exons)) {
    exons <- exon_strip_tx_order(tx_exons[[id]])
    n <- length(exons)
    if (!n) next
    map <- explain_tx_display_map(exons, display)
    if (!nrow(map)) next
    w <- width(exons)
    tx_end <- cumsum(w); tx_start <- tx_end - w + 1L
    p <- explain_map_ranges(tx_start, tx_end, map)
    if (!nrow(p)) next
    lane <- lane + 1L
    k <- p$id
    ex_out[[id]] <- data.table(x0 = p$x0, x1 = p$x1, lane = lane, tx = id, exon = k, n_exons = n,
                               parity = ifelse(k %% 2L == 1L, "odd exon", "even exon"),
                               seqnames = as.character(seqnames(exons))[k], start = start(exons)[k],
                               end = end(exons)[k], strand = as.character(strand(exons))[k],
                               width = w[k], tx_start = tx_start[k], tx_end = tx_end[k],
                               shown = p$x1 - p$x0 + 1L)
    if (isTRUE(collapsed) && n > 1L) {
      # connector between consecutive exons that are both in view and separated on the display
      e <- ex_out[[id]]
      e <- e[, .(x0 = min(x0), x1 = max(x1)), by = exon]
      nxt <- match(e$exon + 1L, e$exon)
      ok <- which(!is.na(nxt))
      if (length(ok)) {
        a <- e[ok]; b <- e[nxt[ok]]
        lo <- pmin(a$x1, b$x1); hi <- pmax(a$x0, b$x0)
        gap <- hi - lo - 1L
        keep <- gap > 0L
        if (any(keep)) {
          ki <- a$exon[keep]
          g0 <- pmin(end(exons)[ki], end(exons)[ki + 1L]) + 1L
          g1 <- pmax(start(exons)[ki], start(exons)[ki + 1L]) - 1L
          in_out[[id]] <- data.table(x0 = lo[keep] + 1L, x1 = hi[keep] - 1L, lane = lane, tx = id,
                                     intron = ki, seqnames = as.character(seqnames(exons))[ki],
                                     start = g0, end = g1, width = g1 - g0 + 1L)
        }
      }
    }
  }
  if (!length(ex_out)) return(NULL)
  list(exons = rbindlist(ex_out), introns = rbindlist(in_out), n_lanes = lane)
}

#' Plotly panel of the exon strip (one thin lane per transcript).
#' @noRd
exon_strip_panel <- function(sd, display_range) {
  dt <- sd$exons
  nl <- sd$n_lanes
  cols <- exon_strip_colors()
  base <- nl - dt$lane + 0.15
  hover <- sprintf("Exon %d of %d (%s)<br>%s:%s-%s (%s)<br>%s nt; transcript nt %s-%s%s",
                   dt$exon, dt$n_exons, dt$tx, dt$seqnames,
                   format(dt$start, big.mark = ",", trim = TRUE), format(dt$end, big.mark = ",", trim = TRUE),
                   dt$strand, format(dt$width, big.mark = ",", trim = TRUE),
                   format(dt$tx_start, big.mark = ",", trim = TRUE), format(dt$tx_end, big.mark = ",", trim = TRUE),
                   ifelse(dt$shown < dt$width, sprintf("<br>(%d nt in view)", dt$shown), ""))
  traces <- explain_bar_traces(dt, dt$parity, cols, rep(0.7, nrow(dt)), base, hover)
  for (k in seq_along(traces)) traces[[k]]$meta <- list(rc_track = "exon_strip")
  it <- sd$introns
  if (!is.null(it) && nrow(it)) {
    y <- nl - it$lane + 0.5
    ih <- sprintf("Intron %d (%s, collapsed)<br>%s:%s-%s<br>%s nt", it$intron, it$tx, it$seqnames,
                  format(it$start, big.mark = ",", trim = TRUE), format(it$end, big.mark = ",", trim = TRUE),
                  format(it$width, big.mark = ",", trim = TRUE))
    traces[[length(traces) + 1L]] <- list(
      type = "scatter", mode = "lines", x = I(as.vector(rbind(it$x0 - 0.5, it$x1 + 0.5, NA))),
      y = I(as.vector(rbind(y, y, NA))), text = I(as.vector(rbind(ih, ih, NA))), hoverinfo = "text",
      line = list(color = "#8a8a8a", width = 1.5), connectgaps = FALSE, showlegend = FALSE,
      name = "collapsed intron", xaxis = "x", yaxis = "y", meta = list(rc_track = "exon_strip_intron"))
  }
  # exon numbers on blocks wide enough to hold them (secondary encoding next to the colour)
  L <- sum(width(if (is(display_range, "GRangesList")) unlist(display_range, use.names = FALSE) else display_range))
  wide <- (dt$x1 - dt$x0 + 1) >= max(1, L / 80)
  if (any(wide)) {
    lab <- dt[wide]
    lcol <- exon_strip_label_colors()
    traces[[length(traces) + 1L]] <- list(
      type = "scatter", mode = "text", x = I((lab$x0 + lab$x1) / 2), y = I(nl - lab$lane + 0.5),
      text = I(as.character(lab$exon)), textfont = list(size = 9, color = I(unname(lcol[lab$parity]))),
      hoverinfo = "skip", showlegend = FALSE, xaxis = "x", yaxis = "y", name = "exon number",
      meta = list(rc_track = "exon_strip_labels"))
  }
  explain_panel(traces, display_range, "Exons", c(0, nl))
}

#' Exon strip controller: which transcripts get a lane (plain list, cache safe).
#'
#' Active when the browser setting is on or explain mode is on with the layer "exons" selected (also
#' for transcripts without explainer results); lanes = the selected transcript, plus the
#' explainability transcripts when explain mode has results at the locus.
#' @noRd
exon_strip_controller <- function(input, tx_grl, explain = NULL, collapsed_introns_width = 0) {
  explain_on <- isTRUE(browser_input_or_default(input, "explain_mode", FALSE)) &&
    "exons" %in% browser_input_or_default(input, "explain_layers", c("status", "expected", "exons"))
  if (!isTRUE(browser_input_or_default(input, "exon_strip", FALSE)) && !explain_on) {
    return(list(active = FALSE))
  }
  ids <- shiny::isolate(input$tx)
  if (explain_on && isTRUE(explain$active)) ids <- c(ids, explain$ids)
  ids <- head(unique(intersect(ids, names(tx_grl))), 5L)
  if (!length(ids)) return(list(active = FALSE))
  list(active = TRUE, ids = ids, tx_exons = tx_grl[ids],
       collapsed = isTRUE(collapsed_introns_width > 0) && isTRUE(browser_input_or_default(input, "viewMode", FALSE)),
       after_status = explain_on)
}

#' Exon strip panel for the controller state, NULL when off / nothing in view.
#' @noRd
exon_strip_track_panel <- function(exon_strip, display_range) {
  if (!isTRUE(exon_strip$active)) return(NULL)
  sd <- exon_strip_data(exon_strip$tx_exons, display_range, exon_strip$collapsed)
  if (is.null(sd)) return(NULL)
  exon_strip_panel(sd, display_range)
}

#' Insert the exon strip panel into the custom bottom panels: directly below the explain status
#' lane when present, else after the sequence tracks (before the other explain lanes).
#' @noRd
exon_strip_insert <- function(panels, strip) {
  if (is.null(strip)) return(panels)
  nms <- names(panels)
  if (is.null(nms)) nms <- rep("", length(panels))
  pos <- match("explain_status", nms)
  if (is.na(pos)) {
    ex <- which(startsWith(nms, "explain"))
    pos <- if (length(ex)) ex[1] - 1L else length(panels)
  }
  append(panels, list(exon_strip = strip), after = pos)
}

#' Colour legend of the exon strip (HTML).
#' @noRd
exon_strip_legend_ui <- function() {
  cols <- exon_strip_colors()
  items <- lapply(names(cols), function(s) tags$span(
    style = "display:inline-block; margin-right:8px; white-space:nowrap;",
    tags$span(style = paste0("display:inline-block; width:12px; height:12px; margin-right:3px; ",
                             "vertical-align:middle; background:", cols[[s]], ";")), s))
  tags$div(class = "exon-strip-legend", style = "font-size:12px;",
           tags$b("Exons (numbered 5'→3' of the transcript): "), items)
}
