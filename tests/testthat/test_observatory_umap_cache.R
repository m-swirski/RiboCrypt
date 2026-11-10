test_that("precomputed UMAP cache is rejected for stale sources or collections", {
  p <- structure(list(),class="plotly")
  x <- list(schema=6L,source_md5="old",plotly_version="1",runs=c("A","B"),plot=p)
  expect_true(rc_umap_cache_valid(x,"old",c("A","B"),"1"))
  expect_false(rc_umap_cache_valid(x,"changed",c("A","B"),"1"))
  expect_false(rc_umap_cache_valid(x,"old",c("B","A"),"1"))
  expect_false(rc_umap_cache_valid(x,"old",c("A","B"),"2"))
  expect_null(rc_umap_cached_plot(NULL,c("tissue","cell_line"),"2d",c(1,2)))
})

test_that("UMAP axes use one projection and do not modify the caller's data", {
  x <- data.table::data.table(`UMAP 1`=1:3,`UMAP 2`=4:6,`UMAP 3`=7:9,
    color_column=c("A","B","A"),BioProject="P",sample=c("A","B","C"),
    author="a",inhibitors="none",tissues_cell_lines="t | c")
  original <- data.table::copy(x)
  p <- umap_plot(x,color.by=c("tissue","cell_line"),template=observatoryUmapPlotlyTemplate(),axes=c(3,2))
  expect_identical(x,original)
  expect_true(all(unlist(lapply(p$x$data, `[[`, "x")) %in% 7:9))
  expect_identical(p$x$layout$xaxis$title,"UMAP 3")
  expect_error(umap_plot(x,axes=c(1,1)),"different")
  expect_error(umap_plot(x,axes=c(1,6)),"not present")
})

test_that("a stale table fallback never exposes samples outside the collection", {
  x <- data.table::data.table(sample=c("a","b","c"),`UMAP 1`=1:3)
  expect_identical(rc_umap_align_runs(x,c("a","b","c")),x)
  expect_warning(y <- rc_umap_align_runs(x,c("b","a","new")),"samples differ")
  expect_identical(y$sample,c("b","a"))
  expect_error(rc_umap_align_runs(x,"not-present"),"no samples")
})

test_that("a missing provenance hash cannot validate a cache", {
  p <- structure(list(), class = "plotly")
  x <- list(schema=6L, source_md5=NA_character_, plotly_version="1", runs="A", plot=p)
  expect_false(rc_umap_cache_valid(x, NA_character_, "A", "1"))
})

test_that("axis annotations follow the selected fit and fail closed when stale", {
  dir <- tempfile(); dir.create(dir); on.exit(unlink(dir, recursive=TRUE))
  source <- file.path(dir, "embedding.fst"); writeLines("coordinates", source)
  annotation <- file.path(dir, "annotations.csv")
  a <- data.frame(axis=1:3, label=c("A +", "B -", "C +"),
                  source_md5=unname(tools::md5sum(source)))
  write.csv(a, annotation, row.names=FALSE)
  Sys.setFileTime(annotation, Sys.time()+2)
  p <- observatoryUmapPlotlyTemplate()
  q <- plotly::plotly_build(rc_umap_apply_axis_annotations(p, source, annotation, 3L, c(2L,3L)))
  expect_identical(q$x$layout$xaxis$title, "UMAP 2 (Enriched: B -)")
  expect_identical(q$x$layout$yaxis$title, "UMAP 3 (Enriched: C +)")
  expect_identical(rc_umap_apply_axis_annotations(p, source, annotation, 4L, c(2L,3L)), p)
  a$source_md5 <- "wrong"; write.csv(a, annotation, row.names=FALSE)
  expect_identical(rc_umap_apply_axis_annotations(p, source, annotation, 3L, c(2L,3L)), p)
  a$source_md5 <- unname(tools::md5sum(source)); write.csv(a, annotation, row.names=FALSE)
  Sys.setFileTime(annotation, Sys.time()-10)
  expect_identical(rc_umap_apply_axis_annotations(p, source, annotation, 3L, c(2L,3L)), p)
})


test_that("prebuilt UMAP widgets skip only the redundant Plotly rebuild", {
  p <- observatoryUmapPlotlyTemplate()
  p <- plotly::plotly_build(p)
  class(p) <- c("rc_umap_prebuilt", class(p))
  expect_identical(plotly::plotly_build(p), p)
  expect_s3_class(p, "plotly")
})

test_that("compact UMAP hover retains the Run and every metadata field", {
  x <- data.table::data.table(`UMAP 1`=1, `UMAP 2`=2,
    color_column="kidney | HEK", BioProject="PRJ123", sample="SRR1",
    author="Smith", inhibitors="chx", tissues_cell_lines="kidney | HEK")
  p <- umap_plot(x, color.by=c("tissue","cell_line"),
                 template=observatoryUmapPlotlyTemplate())
  t <- p$x$data[[1]]
  expect_identical(t$text, "SRR1")
  expect_equal(unname(t$customdata[1,]), c("PRJ123", "Smith", "chx", "kidney | HEK"))
  expect_match(t$hovertemplate, "Run ID: %{text}", fixed=TRUE)
  expect_identical(t$type, "scattergl")
})

test_that("SVG compatibility conversion does not mutate the cached widget", {
  old <- options(ribocrypt.umap_webgl=FALSE); on.exit(options(old), add=TRUE)
  p <- observatoryUmapPlotlyTemplate()
  q <- rc_umap_svg_fallback(p)
  expect_identical(p$x$data[[1]]$type, "scattergl")
  expect_identical(q$x$data[[1]]$type, "scatter")
})
