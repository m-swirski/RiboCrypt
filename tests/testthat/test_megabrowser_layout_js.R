test_that("MegaBrowser layout handles document-targeted Shiny updates", {
  node <- Sys.which("node")
  skip_if(!nzchar(node), "Node.js is required for JavaScript callback tests")
  withr::local_envvar(RIBOCRYPT_MEGA_LAYOUT_JS = system.file("js", "megabrowser_layout.js", package = "RiboCrypt"))
  result <- suppressWarnings(system2(node, c("--test", shQuote(test_path("js", "megabrowser-layout.cjs"))),
    stdout = TRUE, stderr = TRUE))
  expect_null(attr(result, "status"), info = paste(result, collapse = "\n"))
})
