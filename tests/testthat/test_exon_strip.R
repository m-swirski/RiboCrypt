library(GenomicRanges)
library(data.table)

# plus: exons 1001-1100, 1401-1450, 1801-1900 (250 nt); minus: same coordinates on "-"
tx_fixture <- function() {
  GRangesList(
    TXP = GRanges("1", IRanges(c(1001, 1401, 1801), c(1100, 1450, 1900)), "+"),
    TXM = GRanges("1", IRanges(c(1801, 1401, 1001), c(1900, 1450, 1100)), "-"))
}

test_that("exons are numbered 5'->3' on both strands, whatever the input order", {
  grl <- tx_fixture()
  # transcript view, plus strand
  sd <- RiboCrypt:::exon_strip_data(grl["TXP"], grl["TXP"])
  e <- sd$exons
  expect_equal(e$exon, 1:3)
  expect_equal(e$start, c(1001, 1401, 1801))
  expect_equal(c(e$x0, e$x1), c(1, 101, 151, 100, 150, 250))
  # minus strand: exon 1 = highest coordinate, at the left of the display
  sd <- RiboCrypt:::exon_strip_data(grl["TXM"], grl["TXM"])
  e <- sd$exons
  expect_equal(e$exon, 1:3)
  expect_equal(e$start, c(1801, 1401, 1001))
  expect_equal(e$x0, c(1, 101, 151)); expect_equal(e$x1, c(100, 150, 250))
  expect_equal(e$tx_start, c(1, 101, 151))
  # unsorted minus input gives the same numbering
  rev_in <- GRangesList(TXM = rev(grl[["TXM"]]))
  e2 <- RiboCrypt:::exon_strip_data(rev_in, grl["TXM"])$exons
  expect_equal(e2$start, e$start); expect_equal(e2$x0, e$x0)
})

test_that("colours alternate odd / even and are two distinct neutrals", {
  grl <- tx_fixture()
  e <- RiboCrypt:::exon_strip_data(grl["TXP"], grl["TXP"])$exons
  expect_equal(e$parity, c("odd exon", "even exon", "odd exon"))
  cols <- RiboCrypt:::exon_strip_colors()
  expect_length(unique(cols), 2)
  p <- RiboCrypt:::exon_strip_panel(RiboCrypt:::exon_strip_data(grl["TXP"], grl["TXP"]), grl["TXP"])
  bars <- Filter(function(t) t$type == "bar", p$x$data)
  expect_setequal(vapply(bars, function(t) t$marker$color, ""), unname(cols))
  odd <- Filter(function(t) t$name == "odd exon", bars)[[1]]
  expect_equal(as.numeric(odd$x), c(50.5, 200.5))
  expect_true(any(grepl("Exon 2 of 3 \\(TXP\\)<br>1:1,401-1,450 \\(\\+\\)", unlist(lapply(bars, `[[`, "hovertext")))))
  # exon numbers drawn as text
  labs <- Filter(function(t) identical(t$mode, "text"), p$x$data)
  expect_equal(as.character(labs[[1]]$text), c("1", "2", "3"))
})

test_that("genomic view: blocks over exons with empty introns, on both strands", {
  grl <- tx_fixture()
  span_p <- GRangesList(TXP = GRanges("1", IRanges(1001, 1900), "+"))
  e <- RiboCrypt:::exon_strip_data(grl["TXP"], span_p)$exons
  expect_equal(e$x0, c(1, 401, 801)); expect_equal(e$x1, c(100, 450, 900))
  span_m <- GRangesList(TXM = GRanges("1", IRanges(1001, 1900), "-"))
  sd <- RiboCrypt:::exon_strip_data(grl["TXM"], span_m)
  e <- sd$exons
  expect_equal(e$exon, 1:3)
  expect_equal(e$x0, c(1, 451, 801)); expect_equal(e$x1, c(100, 500, 900))
  expect_equal(nrow(sd$introns), 0)  # not collapsed: no connectors
  # partly in view
  part <- GRangesList(TXP = GRanges("1", IRanges(1051, 1420), "+"))
  e <- RiboCrypt:::exon_strip_data(grl["TXP"], part)$exons
  expect_equal(e$exon, 1:2); expect_equal(e$shown, c(50, 20))
})

test_that("collapsed introns: exon blocks plus connectors make the strip contiguous", {
  grl <- tx_fixture()
  for (id in c("TXP", "TXM")) {
    # pseudo-intron display: exons with 30 nt intron flanks
    d <- reduce(c(grl[[id]], flank(grl[[id]], 30, start = TRUE), flank(grl[[id]], 30, start = FALSE)))
    d <- d[order(start(d), decreasing = id == "TXM")]
    sd <- RiboCrypt:::exon_strip_data(grl[id], GRangesList(d), collapsed = TRUE)
    e <- sd$exons; it <- sd$introns
    expect_equal(e$exon, 1:3)
    expect_equal(nrow(it), 2)
    expect_equal(it$intron, 1:2)
    expect_equal(it$width, c(300, 350)[if (id == "TXP") 1:2 else 2:1])
    # exon, connector, exon, connector, exon tile the span without gaps
    segs <- rbind(e[, .(x0, x1)], it[, .(x0, x1)])
    segs <- segs[order(x0)]
    expect_true(all(segs$x0[-1] == segs$x1[-nrow(segs)] + 1L))
    p <- RiboCrypt:::exon_strip_panel(sd, GRangesList(d))
    expect_true(any(vapply(p$x$data, function(t) identical(t$meta$rc_track, "exon_strip_intron"), logical(1))))
  }
})

test_that("controller, insertion below explain status and custom panel polish", {
  grl <- tx_fixture()
  input <- list(tx = "TXM", exon_strip = FALSE, viewMode = TRUE)
  expect_false(RiboCrypt:::exon_strip_controller(input, grl)$active)
  input$exon_strip <- TRUE
  ctl <- RiboCrypt:::exon_strip_controller(input, grl, collapsed_introns_width = 30)
  expect_true(ctl$active); expect_equal(ctl$ids, "TXM"); expect_true(ctl$collapsed)
  # explain mode with the exons layer turns it on without the setting
  input$exon_strip <- FALSE
  input$explain_mode <- TRUE
  input$explain_layers <- c("status", "exons")
  ex <- list(active = TRUE, layers = c("status", "exons"), ids = c("TXM", "TXP"))
  ctl <- RiboCrypt:::exon_strip_controller(input, grl, ex)
  expect_equal(ctl$ids, c("TXM", "TXP"))
  # transcript without explainer results (explain inactive): strip still shown for the selected tx
  ctl <- RiboCrypt:::exon_strip_controller(input, grl, list(active = FALSE))
  expect_true(ctl$active); expect_equal(ctl$ids, "TXM")
  input$explain_layers <- "status"
  expect_false(RiboCrypt:::exon_strip_controller(input, grl, ex)$active)

  panels <- list(phylo = 1, explain_status = 2, explain_plddt = 3)
  expect_equal(names(RiboCrypt:::exon_strip_insert(panels, 9)), c("phylo", "explain_status", "exon_strip", "explain_plddt"))
  expect_equal(names(RiboCrypt:::exon_strip_insert(list(phylo = 1, explain_plddt = 3), 9)),
               c("phylo", "exon_strip", "explain_plddt"))
  expect_equal(names(RiboCrypt:::exon_strip_insert(list(), 9)), "exon_strip")
  expect_equal(RiboCrypt:::exon_strip_insert(panels, NULL), panels)

  panel <- RiboCrypt:::exon_strip_track_panel(ctl, grl["TXM"])
  out <- RiboCrypt:::custom_panels_to_plotly(list(exon_strip = panel))
  expect_identical(out$exon_strip, panel)
  expect_true("exons" %in% RiboCrypt:::explain_layer_choices())
})
