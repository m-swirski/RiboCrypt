devtools::load_all(".")
root <- Sys.getenv("RIBOCRYPT_REPORT_TEST_ROOT")
stopifnot(nzchar(root), dir.exists(root))
database <- file.path(root, "reports.sqlite")
unlink(database)
con <- ribocrypt_access_db(database)
ribocrypt_access_workspace(con, "alice")
ribocrypt_access_member(con, "issuer", "alice", "alice")
ribocrypt_access_dataset(con, "private", "private", root)
ribocrypt_access_grant(con, "alice", "private")
DBI::dbDisconnect(con)
writeLines("<html><body>SELECTED_REPORT<script>document.body.dataset.script='works'</script></body></html>",
           file.path(root, "sample.html"))
writeLines("PRIVATE_SIBLING", file.path(root, "sibling.txt"))
config <- ribocrypt_access_control(database, "issuer", strrep("x", 64), refresh_ms = 60000)
app <- shiny::shinyApp(shiny::fluidPage(shiny::uiOutput("report")), function(input, output, session) {
  # Chromium does not send extraHTTPHeaders on WebSocket handshakes. This fixed
  # synthetic cookie is ONLY for this disposable loopback fixture, never production.
  request <- session$request
  if (grepl("rc_report_test=alice", request$HTTP_COOKIE %||% "", fixed = TRUE)) {
    request$HTTP_RIBOCRYPT_GATEWAY_SECRET <- strrep("x", 64)
    request$HTTP_RIBOCRYPT_SUBJECT <- "alice"
  }
  context <- RiboCrypt:::rc_access_attach(session, config, request)
  if (!"private" %in% context$catalog$id) return()
  # Exercise the real selection helper, mocking only the ORFik experiment fixture.
  url <- testthat::with_mocked_bindings(
    RiboCrypt:::get_fastq_page(list(library = "sample"), NULL, NULL, "."),
    observed_exp_subset = function(...) list(filepath = file.path(root, "sample.bam")),
    libFolder = function(...) root, .package = "RiboCrypt")
  output$report <- shiny::renderUI(shiny::tags$iframe(id = "report-frame", src = url,
                                                    sandbox = "allow-scripts"))
})
shiny::runApp(app, host = "127.0.0.1", port = 7849, launch.browser = FALSE)
