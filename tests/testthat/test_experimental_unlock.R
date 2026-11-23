# Session unlock for experimental / private features (dummy codes only).

dummy <- "test-unlock-value"

reset_global_attempts <- function() {
  assign("wrong", numeric(), envir = RiboCrypt:::.rc_experimental_env)
}

clean_unlock_config <- function(env = parent.frame()) {
  withr::local_options(RiboCrypt.experimental_unlock = NULL,
                       RiboCrypt.experimental_unlock_file = NULL, .local_envir = env)
  withr::local_envvar(RIBOCRYPT_EXPERIMENTAL_UNLOCK = NA,
                      RIBOCRYPT_EXPERIMENTAL_UNLOCK_FILE = NA, .local_envir = env)
}

has_bytes <- function(x, pattern) {
  length(grepRaw(pattern, x, fixed = TRUE)) > 0
}

# Module that mimics the browser / observatory module wiring.
unlock_module <- function(id, verifier, clock) {
  shiny::moduleServer(id, function(input, output, session) {
    RiboCrypt:::experimental_session_init(session, verifier, now = function() clock$t)
    RiboCrypt:::experimental_unlock_server(input, output, session)
  })
}

test_that("locked by default: no configuration means no verifier and no gate", {
  clean_unlock_config()
  expect_null(RiboCrypt:::experimental_unlock_setup(NULL))
  expect_false(RiboCrypt:::experimental_feature_enabled(NULL, "explain"))
  expect_null(RiboCrypt:::experimental_unlock_ui(identity, configured = FALSE))
  expect_equal(RiboCrypt:::experimental_features_key(NULL), "")
})

test_that("verifier accepts only the right code and holds no plain copy", {
  clean_unlock_config()
  v <- RiboCrypt:::experimental_unlock_setup(dummy)
  expect_s3_class(v, "rc_experimental_verifier")
  expect_true(v$verify(dummy))
  expect_false(v$verify("wrong-value"))
  expect_false(v$verify(paste0(dummy, " ")))
  expect_false(v$verify(""))
  expect_false(v$verify(NA_character_))
  expect_false(v$verify(NULL))
  expect_false(v$verify(c(dummy, dummy)))
  expect_false(v$verify(strrep("x", 5000)))
  expect_false(has_bytes(serialize(v, NULL), dummy))
  expect_false(grepl(dummy, paste(capture.output(print(v)), collapse = "\n"), fixed = TRUE))
  # passing a verifier again keeps it (access-control rebuilds the app per session)
  expect_identical(RiboCrypt:::experimental_unlock_setup(v), v)
})

test_that("option and environment values are read, then cleared; errors never echo the value", {
  clean_unlock_config()
  options(RiboCrypt.experimental_unlock = dummy)
  v <- RiboCrypt:::experimental_unlock_setup()
  expect_true(v$verify(dummy))
  expect_null(getOption("RiboCrypt.experimental_unlock"))

  Sys.setenv(RIBOCRYPT_EXPERIMENTAL_UNLOCK = dummy)
  v <- RiboCrypt:::experimental_unlock_setup()
  expect_true(v$verify(dummy))
  expect_equal(Sys.getenv("RIBOCRYPT_EXPERIMENTAL_UNLOCK"), "")

  err <- tryCatch(RiboCrypt:::experimental_unlock_setup(c(dummy, "b")), error = conditionMessage)
  expect_false(grepl(dummy, err, fixed = TRUE))
})

test_that("unlock code from a file (option, env var and argument)", {
  clean_unlock_config()
  f <- tempfile("rc_unlock")
  writeLines(dummy, f)
  Sys.chmod(f, "0600")
  expect_silent(v <- RiboCrypt:::experimental_unlock_setup(RiboCrypt::experimental_unlock_file(f)))
  expect_true(v$verify(dummy))
  options(RiboCrypt.experimental_unlock_file = f)
  expect_true(RiboCrypt:::experimental_unlock_setup()$verify(dummy))
  options(RiboCrypt.experimental_unlock_file = NULL)
  Sys.setenv(RIBOCRYPT_EXPERIMENTAL_UNLOCK_FILE = f)
  expect_true(RiboCrypt:::experimental_unlock_setup()$verify(dummy))
  Sys.unsetenv("RIBOCRYPT_EXPERIMENTAL_UNLOCK_FILE")
  Sys.chmod(f, "0644")
  expect_warning(RiboCrypt:::experimental_unlock_setup(RiboCrypt::experimental_unlock_file(f)),
                 "chmod 600")
  expect_false(grepl(dummy, paste(capture.output(print(RiboCrypt::experimental_unlock_file(f))),
                                  collapse = ""), fixed = TRUE))
  expect_error(RiboCrypt:::experimental_unlock_setup(RiboCrypt::experimental_unlock_file(tempfile())),
               "does not exist")
})

test_that("testServer: locked, wrong code rejected, right code unlocks the session", {
  clean_unlock_config()
  reset_global_attempts()
  clock <- new.env(); clock$t <- 1000
  v <- RiboCrypt:::experimental_unlock_setup(dummy)
  shiny::testServer(unlock_module, args = list(verifier = v, clock = clock), {
    expect_false(RiboCrypt:::experimental_feature_enabled(session, "explain"))
    expect_error(output$explain_settings)  # req() fails: nothing rendered while locked

    session$setInputs(experimental_unlock_code = "wrong-value", experimental_unlock_go = 1)
    expect_equal(output$experimental_unlock_status, "Wrong code")
    expect_false(RiboCrypt:::experimental_feature_enabled(session, "explain"))

    clock$t <- clock$t + 5
    session$setInputs(experimental_unlock_code = dummy, experimental_unlock_go = 2)
    expect_equal(output$experimental_unlock_status, "Unlocked")
    expect_true(RiboCrypt:::experimental_feature_enabled(session, "explain"))
    expect_equal(RiboCrypt:::experimental_features_key(session), "explain")
    # outputs sent to the client never contain the code
    expect_false(grepl(dummy, output$experimental_unlock_status, fixed = TRUE))
    expect_false(grepl(dummy, paste(unlist(output$explain_settings), collapse = ""), fixed = TRUE))
    # session state holds no plain copy
    expect_false(has_bytes(serialize(session$userData$rc_experimental$verify, NULL), dummy))
  })
  # a new session starts locked again
  shiny::testServer(unlock_module, args = list(verifier = v, clock = clock), {
    expect_false(RiboCrypt:::experimental_feature_enabled(session, "explain"))
  })
})

test_that("no unlock configured: typing any code never unlocks", {
  reset_global_attempts()
  clock <- new.env(); clock$t <- 1000
  shiny::testServer(unlock_module, args = list(verifier = NULL, clock = clock), {
    session$flushReact()
    session$setInputs(experimental_unlock_code = dummy, experimental_unlock_go = 1)
    expect_equal(output$experimental_unlock_status, "Wrong code")
    expect_false(RiboCrypt:::experimental_feature_enabled(session, "explain"))
  })
})

test_that("rate limit: cool-down after a wrong code and max attempts per minute", {
  clean_unlock_config()
  reset_global_attempts()
  clock <- new.env(); clock$t <- 1000
  v <- RiboCrypt:::experimental_unlock_setup(dummy)
  shiny::testServer(unlock_module, args = list(verifier = v, clock = clock), {
    session$flushReact()
    try_code <- function(code, n) {
      session$setInputs(experimental_unlock_code = code, experimental_unlock_go = n)
      output$experimental_unlock_status
    }
    expect_equal(try_code("bad", 1), "Wrong code")
    clock$t <- clock$t + 1  # within the cool-down: even the right code waits
    expect_match(try_code(dummy, 2), "Too many attempts")
    expect_false(RiboCrypt:::experimental_feature_enabled(session, "explain"))
    for (i in 3:6) { clock$t <- clock$t + 3; expect_equal(try_code("bad", i), "Wrong code") }
    clock$t <- clock$t + 3  # 5 wrong codes in the last minute
    expect_match(try_code(dummy, 7), "Too many attempts")
    expect_false(RiboCrypt:::experimental_feature_enabled(session, "explain"))
    clock$t <- clock$t + 61
    expect_equal(try_code(dummy, 8), "Unlocked")
  })
})

test_that("process-wide limit stops bots that open many sessions", {
  clean_unlock_config()
  reset_global_attempts()
  v <- RiboCrypt:::experimental_unlock_setup(dummy)
  t0 <- 5000
  for (i in seq_len(50)) {
    s <- shiny::MockShinySession$new()
    RiboCrypt:::experimental_session_init(s, v, now = function() t0 + i / 100)
    expect_equal(RiboCrypt:::experimental_try_unlock(s, "bad"), "wrong")
  }
  s <- shiny::MockShinySession$new()
  RiboCrypt:::experimental_session_init(s, v, now = function() t0 + 1)
  expect_equal(RiboCrypt:::experimental_try_unlock(s, dummy), "wait")
  s2 <- shiny::MockShinySession$new()
  RiboCrypt:::experimental_session_init(s2, v, now = function() t0 + 70)
  expect_equal(RiboCrypt:::experimental_try_unlock(s2, dummy), "unlocked")
  reset_global_attempts()
})

test_that("UI: masked password box only when configured, never carries a value", {
  ui <- RiboCrypt:::experimental_unlock_ui(shiny::NS("browser"), configured = TRUE)
  html <- as.character(ui)
  expect_match(html, 'type="password"', fixed = TRUE)
  expect_match(html, 'id="browser-experimental_unlock_code"', fixed = TRUE)
  expect_false(grepl(dummy, html, fixed = TRUE))
  expect_null(RiboCrypt:::experimental_unlock_ui(shiny::NS("browser"), configured = FALSE))
  # settings: explain settings are a server-rendered output, not static UI
  html2 <- as.character(RiboCrypt:::experimental_settings_ui(shiny::NS("browser")))
  expect_match(html2, 'id="browser-explain_settings"', fixed = TRUE)
  expect_false(grepl("Explainability", html2, fixed = TRUE))
})

test_that("share URL and cache key never contain the code; explain gated by unlock", {
  input <- list(dff = "exp", gene = "G", tx = "TX1", explain_mode = TRUE,
                experimental_unlock_code = dummy, experimental_unlock_go = 1)
  url <- RiboCrypt:::make_url_from_inputs_parameters(input)
  expect_false(grepl(dummy, url, fixed = TRUE))
  expect_false(grepl("experimental", url, fixed = TRUE))
  expect_false(grepl("explain_mode", url, fixed = TRUE))  # locked: not shared

  rv <- shiny::reactiveValues(gene = "G", tx = "TX1", experimental_unlock_code = dummy,
                              experimental_unlock_go = 1)
  key <- shiny::isolate(RiboCrypt:::input_to_list(rv, list(userAgent = "x")))
  expect_false(any(grepl("experimental_unlock", names(key), fixed = TRUE)))
  expect_false(grepl(dummy, paste(unlist(key), collapse = ""), fixed = TRUE))
  expect_equal(key$rc_experimental, "")

  # explain hash and controller depend on the unlock state
  inp <- list(explain_mode = TRUE, explain_layers = "status")
  expect_false(identical(RiboCrypt:::explain_hash(inp, enabled = FALSE),
                         RiboCrypt:::explain_hash(inp, enabled = TRUE)))
  expect_false(RiboCrypt:::explain_controller(inp, NULL, NULL, NULL)$active)
})

test_that("observatory share URL settings never include the code", {
  input <- list(experimental_unlock_code = dummy, gene = "G", tx = "TX1")
  settings <- RiboCrypt:::observatory_capture_browser_settings(input)
  expect_false(grepl(dummy, paste(unlist(settings), collapse = ""), fixed = TRUE))
})
