test_that("smoothing extensions are shared only for repeated k-mers", {
  gr <- GenomicRanges::GRangesList(tx = GenomicRanges::GRanges(
    "chr1", IRanges::IRanges(c(1000, 1200), c(1099, 1299)), strand = "+"
  ))
  ranges <- RiboCrypt:::browser_shared_smoothing_ranges(gr, c(1, 3, 3, 5), 4)
  expect_null(ranges[[1]])
  expect_null(ranges[[4]])
  expected <- ORFik::extendTrailers(ORFik::extendLeaders(gr, 9), 9)
  expect_identical(ranges[[2]], expected)
  expect_identical(ranges[[3]], expected)
  expect_identical(RiboCrypt:::browser_shared_smoothing_ranges(gr, 1, 2),
                   list(NULL, NULL))
  expect_identical(RiboCrypt:::browser_shared_smoothing_ranges(gr, 3, 1),
                   list(NULL))
})

test_that("shared smoothing preserves each coverage profile", {
  df <- ORFik::ORFik.template.experiment()[9, ]
  gr <- ORFik::loadRegion(df, "mrna")[1]
  reads <- ORFik::filepath(df, "bigwig")
  expected <- list(
    RiboCrypt:::getProfileWrapper(gr, reads[[1]], TRUE, 3, FALSE),
    RiboCrypt:::getProfileWrapper(gr, reads[[1]], FALSE, 3, TRUE)
  )
  got <- RiboCrypt:::build_browser_profiles(
    gr, rep(reads, 2), kmers = 3, kmers_type = "mean", frames_type = "area",
    frames_subset = "all", withFrames = c(TRUE, FALSE),
    log_scale = c(FALSE, TRUE), BPPARAM = BiocParallel::SerialParam()
  )
  expect_identical(got, expected)
  if (.Platform$OS.type != "windows") {
    parallel <- RiboCrypt:::build_browser_profiles(
      gr, rep(reads, 2), kmers = 3, kmers_type = "mean", frames_type = "area",
      frames_subset = "all", withFrames = c(TRUE, FALSE),
      log_scale = c(FALSE, TRUE), BPPARAM = BiocParallel::MulticoreParam(2)
    )
    expect_equal(parallel, expected)
    for (i in seq_along(expected)) {
      expect_identical(names(parallel[[i]]), names(expected[[i]]))
      for (column in names(expected[[i]])) {
        expect_identical(parallel[[i]][[column]], expected[[i]][[column]])
      }
    }
  }
})
