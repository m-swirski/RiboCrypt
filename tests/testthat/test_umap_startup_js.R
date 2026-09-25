test_that("the shipped UMAP callback skips no-op updates without losing selection resets", {
  node <- Sys.which("node")
  skip_if(!nzchar(node), "Node.js is required for JavaScript callback tests")
  withr::local_envvar(RIBOCRYPT_UMAP_JS = system.file("js", "umap_plot_extension.js", package = "RiboCrypt"))
  result <- suppressWarnings(system2(node, c("--test", shQuote(test_path("js", "umap-startup.cjs"))),
                                    stdout = TRUE, stderr = TRUE))
  expect_null(attr(result, "status"), info = paste(result, collapse = "\n"))
})
