session_access_fixture <- function() {
  skip_if_not_installed("DBI")
  skip_if_not_installed("RSQLite")
  directory <- tempfile("account-experiments-")
  dir.create(directory)
  template <- ORFik::ORFik.template.experiment(as.temp = TRUE)
  template[[ncol(template) + 1L]] <- ""
  template[4, ncol(template)] <- "Run"
  for (name in c("public-exp", "alice-exp", "bob-exp")) {
    template[1, 2] <- name
    template[5:nrow(template), ncol(template)] <- paste0(name, "-", seq_len(nrow(template) - 4L))
    ORFik::save.experiment(template, file.path(directory, paste0(name, ".csv")))
  }
  database <- file.path(directory, "access.sqlite")
  con <- ribocrypt_access_db(database)
  for (user in c("alice", "bob")) {
    ribocrypt_access_workspace(con, user)
    ribocrypt_access_member(con, "issuer", user, user)
  }
  for (name in c("public-exp", "alice-exp", "bob-exp"))
    ribocrypt_access_dataset(con, name, name, directory, public = name == "public-exp")
  ribocrypt_access_grant(con, "alice", "alice-exp")
  ribocrypt_access_grant(con, "bob", "bob-exp", download = TRUE)
  DBI::dbDisconnect(con)
  list(directory = directory, database = database,
       config = ribocrypt_access_control(database, "issuer", paste(rep("x", 64), collapse = "")))
}

access_session_request <- function(config, user = NULL, path = "/") {
  list(HTTP_RIBOCRYPT_GATEWAY_SECRET = config$gateway_secret,
       HTTP_RIBOCRYPT_SUBJECT = user, PATH_INFO = path, REQUEST_METHOD = "GET", QUERY_STRING = "")
}

test_that("FASTQ reports publish only a session object with live read permission", {
  f <- session_access_fixture()
  on.exit(unlink(f$directory, recursive = TRUE))
  con <- ribocrypt_access_db(f$database)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  identity <- list(issuer = "issuer", subject = "alice")
  ctx <- list(con = con, identity = identity, catalog = ribocrypt_access_catalog(con, identity))
  report <- file.path(f$directory, "sample.html")
  writeLines("<html>selected report</html>", report)
  response <- rc_fastq_report_response(report, f$directory, ctx)
  expect_equal(response$status, 200L)
  expect_match(rawToChar(response$body), "selected report")
  expect_identical(response$headers[["Cache-Control"]], "no-store")
  expect_identical(response$headers[["Content-Security-Policy"]], "sandbox allow-scripts")
  outside <- tempfile(fileext = ".html")
  writeLines("private sibling", outside)
  on.exit(unlink(outside), add = TRUE)
  link <- file.path(f$directory, "escape.html")
  expect_true(file.symlink(outside, link))
  expect_equal(rc_fastq_report_response(link, f$directory, ctx)$status, 403L)
  published <- NULL
  session <- list(userData = list(ribocrypt_access = ctx), registerDataObj = function(name, data, filterFunc) {
    published <<- filterFunc
    "/session/test/dataobj/fastq_report"
  })
  resources <- shiny::resourcePaths()
  expect_match(rc_fastq_report_url(report, f$directory, session), "/dataobj/")
  expect_identical(shiny::resourcePaths(), resources)
  expect_match(rawToChar(published(NULL, list(QUERY_STRING = "file=escape.html"))$body), "selected report")
  ribocrypt_access_grant(con, "alice", "alice-exp", enabled = FALSE)
  expect_equal(rc_fastq_report_response(report, f$directory, ctx)$status, 403L)
  expect_equal(published(NULL, list())$status, 403L)
})

test_that("authenticated app sends only a public shell before identity is checked", {
  f <- session_access_fixture()
  on.exit(unlink(f$directory, recursive = TRUE))
  app <- RiboCrypt_app(access_control = f$config)
  response <- app$httpHandler(access_session_request(f$config))
  expect_equal(response$status, 200)
  expect_false(grepl("alice-exp|bob-exp|public-exp", response$content %||% response$body))
  expect_equal(app$httpHandler(list(PATH_INFO = "/", REQUEST_METHOD = "GET"))$status, 403)
})

test_that("a verified session gets its catalog, filtered metadata and isolated cache", {
  f <- session_access_fixture()
  on.exit(unlink(f$directory, recursive = TRUE))
  results <- list()
  for (user in c("alice", "bob")) {
    mock <- shiny::MockShinySession$new()
    shiny::testServer(function(input, output, session) {
      context <- rc_access_attach(session, f$config, access_session_request(f$config, user))
      prepared <- rc_access_experiment_catalog(context)
      session$userData$ribocrypt_access <- prepared$context
      output$experiment <- shiny::renderText(paste(prepared$all_exp$name, collapse = ","))
      cached <- shiny::reactive(user) |> shiny::bindCache("same-key")
      output$value <- shiny::renderText(cached())
    }, {
      expect_false(grepl(if (user == "alice") "bob-exp" else "alice-exp", output$experiment))
      expect_identical(output$value, user)
      expect_error(rc_read_experiment(if (user == "alice") "bob-exp" else "alice-exp"), "access denied")
      expect_error(rc_read_experiment(file.path(f$directory, "public-exp.csv")), "access denied")
      ctx <- session$userData$ribocrypt_access
      expect_false(identical(ORFik::envExp(rc_read_experiment("public-exp")), .GlobalEnv))
      meta <- data.table::data.table(Run = c(ctx$runs[1], "unassigned-run"), TISSUE = c("OK", "secret"))
      expect_identical(rc_access_metadata(meta, ctx)$TISSUE, "OK")
      expect_error(rc_access_reference(rc_read_experiment("public-exp")), "not enabled")
      results[[user]] <<- ctx$catalog$id
    }, session = mock)
  }
  con <- ribocrypt_access_db(f$database)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  expect_setequal(DBI::dbReadTable(con, "accounts")$subject, c("alice", "bob"))
  expect_setequal(results$alice, c("public-exp", "alice-exp"))
  expect_setequal(results$bob, c("public-exp", "bob-exp"))
})

test_that("session HTTP requests reject other identities and read-only downloads", {
  f <- session_access_fixture()
  on.exit(unlink(f$directory, recursive = TRUE))
  con <- ribocrypt_access_db(f$database)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  identity <- list(issuer = "issuer", subject = "alice")
  ctx <- list(con = con, identity = identity, catalog = ribocrypt_access_catalog(con, identity))
  calls <- 0L
  handler <- rc_access_request_handler(function(req) { calls <<- calls + 1L; list(status = 200L) }, ctx, f$config)
  expect_equal(handler(access_session_request(f$config, "alice", "/dataobj/table"))$status, 200)
  expect_equal(handler(access_session_request(f$config, "bob", "/dataobj/table"))$status, 403)
  expect_equal(handler(access_session_request(f$config, "alice", "/download/coverage"))$status, 403)
  expect_equal(handler(access_session_request(f$config, NULL, "/dataobj/table"))$status, 403)
  ribocrypt_access_grant(con, "alice", "alice-exp", enabled = FALSE)
  expect_equal(handler(access_session_request(f$config, "alice", "/dataobj/table"))$status, 403)
  expect_equal(calls, 1)
})

test_that("the real-session HTTP method can be guarded without leaving it unlocked", {
  f <- session_access_fixture()
  on.exit(unlink(f$directory, recursive = TRUE))
  con <- ribocrypt_access_db(f$database)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  ctx <- list(con = con, identity = NULL, catalog = ribocrypt_access_catalog(con))
  session <- new.env()
  session$handleRequest <- function(req) list(status = 200L)
  lockBinding("handleRequest", session)
  rc_access_protect_session_http(session, ctx, f$config)
  expect_true(bindingIsLocked("handleRequest", session))
  expect_equal(session$handleRequest(access_session_request(f$config))$status, 200)
  expect_equal(session$handleRequest(access_session_request(f$config, "alice"))$status, 403)
})

test_that("reference-wide tables and files cannot escape the session manifest", {
  f <- session_access_fixture()
  on.exit(unlink(f$directory, recursive = TRUE))
  mock <- shiny::MockShinySession$new()
  shiny::testServer(function(input, output, session) {
    rc_access_attach(session, f$config, access_session_request(f$config, "alice"))
  }, {
    ctx <- session$userData$ribocrypt_access
    ctx$runs <- "allowed"
    session$userData$ribocrypt_access <- ctx
    table <- data.table::data.table(sample = c("allowed", "hidden"), value = c(1, 99))
    df <- rc_read_experiment("public-exp")
    testthat::local_mocked_bindings(runIDs = function(df) "allowed", .package = "ORFik")
    expect_identical(rc_access_umap(table, df)$value, 1)
    expect_error(rc_access_umap(data.table::data.table(value = 1), df), "identifiers")
    expect_identical(names(rc_access_collection_columns(data.table::data.table(allowed = 1, hidden = 99))), "allowed")
    expect_error(rc_access_path(file.path(f$directory, "access.sqlite"), file.path(f$directory, "other")), "access denied")
    cache <- rc_access_cache(session$cache, ctx)
    cache$set("value", "private")
    expect_identical(cache$get("value"), "private")
    ribocrypt_access_grant(ctx$con, "alice", "alice-exp", enabled = FALSE)
    expect_error(cache$get("value"), "permissions changed")
  }, session = mock)
})

test_that("account configuration rejects unsafe redirects and malformed polling", {
  f <- session_access_fixture()
  on.exit(unlink(f$directory, recursive = TRUE))
  expect_error(ribocrypt_access_control(f$database, "issuer", f$config$gateway_secret, login_url = "//evil.example/"), "same-site")
  expect_error(ribocrypt_access_control(f$database, "issuer", f$config$gateway_secret, refresh_ms = NA_real_), "finite")
  expect_false(grepl(f$config$gateway_secret, paste(capture.output(print(f$config)), collapse = "\n")))
  expect_identical(rc_access_defaults(c(default_experiment = "hidden", default_gene = "SECRET"), "public-exp", character()),
                   c(default_experiment = "public-exp"))
})

test_that("collection statistics and summary exclude libraries outside the experiment", {
  df <- ORFik::ORFik.template.experiment()[9:10, ]
  df$Run <- c("visible-1", "visible-2")
  path <- tempfile(fileext = ".fst")
  on.exit(unlink(path))
  table <- data.table::data.table(`visible-1` = c(1, 2, 3), `visible-2` = c(3, 5, 7), hidden = c(1e6, 1e6, 1e6))
  fst::write_fst(table, path)
  metadata <- data.table::data.table(Run = c("visible-1", "visible-2", "hidden"),
    TISSUE = c("A", "B", "secret"), BioProject = "study")
  result <- compute_collection_table(path, lib_sizes = NULL, df = df, metadata_field = "TISSUE",
    normalization = "maxNormalized", kmer = 1, metadata = metadata, as_list = TRUE)
  expect_identical(colnames(result$table), c("visible-1", "visible-2"))
  expect_equal(attr(result$table, "summary_cov")$count, c(4, 7, 10))
  expect_length(attr(result$table, "km")$cluster, 2)
  expect_false("hidden" %in% attr(result$metadata_field, "runIDs")$Run)
})
