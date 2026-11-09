test_that("precomputed UMAP cache is rejected for stale sources or collections", {
  p <- structure(list(),class="plotly")
  x <- list(schema=1L,source_md5="old",plotly_version="1",runs=c("A","B"),plot=p)
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
