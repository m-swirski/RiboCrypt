test_that("the shipped sequence callback avoids redundant updates and preserves interactions", {
  node <- Sys.which("node")
  skip_if(!nzchar(node), "Node.js is required for JavaScript callback tests")
  withr::local_envvar(RIBOCRYPT_SEQUENCE_JS = system.file("js", "render_on_zoom.js", package = "RiboCrypt"))
  result <- suppressWarnings(system2(node, c("--test", shQuote(test_path("js", "sequence-startup.cjs"))),
                                    stdout = TRUE, stderr = TRUE))
  expect_null(attr(result, "status"), info = paste(result, collapse = "\n"))
})
