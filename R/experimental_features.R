# Session unlock for private / experimental features.
#
# No user accounts: the app is given one shared unlock code at init (argument
# `experimental_unlock` of RiboCrypt_app(), an option, an environment variable or a
# file). Collaborators type the code in the browser (masked password box in the
# browser / observatory settings) to unlock the experimental features for their
# own Shiny session only.
#
# Safety design:
#  - The plain code is read once at app creation, turned into an HMAC-SHA256 with a
#    random per-process salt and dropped; only salt + digest live in a closure
#    (the verifier). Options / environment variables holding the plain value are
#    cleared after reading.
#  - The typed value travels client -> server only (websocket, TLS on the public
#    app), is hashed with the same salt and compared in constant time; the input
#    box is cleared afterwards and the value is never sent back, logged, cached,
#    put in a URL or in an error message (all failures are reported as
#    "Wrong code").
#  - Unlock state is a session-scoped reactiveVal in session$userData; it is not
#    stored in URLs, bookmarks, cookies or local storage. A reload locks again.
#  - Rate limit: per session at most `max_attempts` wrong codes per minute and a
#    short cool-down after each wrong code; a process-wide limit guards against
#    bots that open many sessions. Nothing blocks the R process (no Sys.sleep).
#  - Plot caches are app-wide, so the unlock state is part of every cache key
#    that can contain experimental layers (input_to_list, explain_hash).

#' Registered experimental features (unlocked together by the code).
#' @noRd
rc_experimental_features <- function() c("explain")

#' Experimental unlock code from a file
#'
#' Use in the launcher instead of a literal value, so the code is neither in the
#' launcher script nor in the process command line. The file should be readable
#' only by the app user (chmod 600) and lie outside any web-served directory.
#' The first line is used (trailing whitespace removed).
#' @param path character, path to the file.
#' @return an object for the `experimental_unlock` argument of
#'  \code{\link{RiboCrypt_app}}.
#' @export
#' @examples
#' # RiboCrypt_app(experimental_unlock = experimental_unlock_file("/run/secrets/rc_unlock"))
experimental_unlock_file <- function(path) {
  stopifnot(is.character(path), length(path) == 1L, !is.na(path), nzchar(path))
  structure(list(path = path), class = "rc_experimental_unlock_file")
}

#' @export
print.rc_experimental_unlock_file <- function(x, ...) {
  cat("RiboCrypt experimental unlock file:", x$path, "\n")
  invisible(x)
}

#' @export
print.rc_experimental_verifier <- function(x, ...) {
  cat("RiboCrypt experimental unlock (configured)\n")
  invisible(x)
}

#' Read the code from a file (first line), warn on loose permissions.
#' @noRd
experimental_unlock_read_file <- function(path) {
  path <- path.expand(path)
  if (!file.exists(path)) stop("experimental unlock file does not exist: ", path, call. = FALSE)
  mode <- file.info(path)$mode
  if (!is.na(mode) && bitwAnd(as.integer(mode), strtoi("077", 8L)) != 0L) {
    warning("experimental unlock file is readable by group/others; use chmod 600: ", path,
            call. = FALSE)
  }
  value <- tryCatch(readLines(path, n = 1L, warn = FALSE),
                    error = function(e) stop("could not read experimental unlock file: ", path,
                                             call. = FALSE))
  if (!length(value)) value <- ""
  sub("[[:space:]]+$", "", value)
}

#' Resolve the unlock configuration into a verifier (or NULL = locked, no box).
#'
#' Precedence: argument > option RiboCrypt.experimental_unlock_file >
#' option RiboCrypt.experimental_unlock > env RIBOCRYPT_EXPERIMENTAL_UNLOCK_FILE >
#' env RIBOCRYPT_EXPERIMENTAL_UNLOCK. Plain-value options / env vars are cleared
#' after reading. Error messages never contain the value.
#' @noRd
experimental_unlock_setup <- function(experimental_unlock = NULL) {
  if (inherits(experimental_unlock, "rc_experimental_verifier")) return(experimental_unlock)
  value <- NULL
  if (inherits(experimental_unlock, "rc_experimental_unlock_file")) {
    value <- experimental_unlock_read_file(experimental_unlock$path)
  } else if (!is.null(experimental_unlock)) {
    value <- experimental_unlock
  } else if (!is.null(f <- getOption("RiboCrypt.experimental_unlock_file", NULL))) {
    value <- experimental_unlock_read_file(f)
  } else if (!is.null(v <- getOption("RiboCrypt.experimental_unlock", NULL))) {
    value <- v
  } else if (nzchar(f <- Sys.getenv("RIBOCRYPT_EXPERIMENTAL_UNLOCK_FILE", ""))) {
    value <- experimental_unlock_read_file(f)
  } else if (nzchar(v <- Sys.getenv("RIBOCRYPT_EXPERIMENTAL_UNLOCK", ""))) {
    value <- v
  }
  # Drop plain copies held in options / environment (the file path options stay).
  if (!is.null(getOption("RiboCrypt.experimental_unlock", NULL))) options(RiboCrypt.experimental_unlock = NULL)
  if (nzchar(Sys.getenv("RIBOCRYPT_EXPERIMENTAL_UNLOCK", ""))) Sys.unsetenv("RIBOCRYPT_EXPERIMENTAL_UNLOCK")
  if (is.null(value)) return(NULL)
  if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(value)) {
    stop("experimental_unlock must be a single non-empty string or experimental_unlock_file()",
         call. = FALSE)
  }
  verifier <- experimental_verifier(value)
  rm(value)
  verifier
}

#' Build a verifier closure holding only a random salt and HMAC-SHA256 digest.
#' @noRd
experimental_verifier <- function(value) {
  salt <- openssl::rand_bytes(32L)
  mac <- function(x) as.raw(openssl::sha256(charToRaw(enc2utf8(x)), key = salt))
  digest <- mac(value)
  rm(value)
  verify <- function(candidate) {
    ok <- tryCatch({
      if (!is.character(candidate) || length(candidate) != 1L || is.na(candidate) ||
          !nzchar(candidate) || nchar(candidate, type = "bytes") > 1024L) {
        FALSE
      } else {
        h <- mac(candidate)
        # constant-time comparison of equal-length digests
        sum(bitwXor(as.integer(h), as.integer(digest))) == 0L
      }
    }, error = function(e) FALSE)
    isTRUE(ok)
  }
  structure(list(verify = verify), class = "rc_experimental_verifier")
}

# Process-wide wrong-attempt log (bots opening many sessions).
.rc_experimental_env <- new.env(parent = emptyenv())

#' Initialise the per-session experimental state (root session).
#' @param session Shiny session.
#' @param verifier result of experimental_unlock_setup() (NULL = not configured).
#' @noRd
experimental_session_init <- function(session, verifier = NULL,
                                      now = function() as.numeric(Sys.time()),
                                      max_attempts = 5L, cooldown = 2,
                                      global_max_attempts = 50L) {
  if (!is.null(session$userData$rc_experimental)) return(invisible(session$userData$rc_experimental))
  state <- list(
    configured = !is.null(verifier),
    unlocked = shiny::reactiveVal(character()),
    verify = if (is.null(verifier)) function(x) FALSE else verifier$verify,
    attempts = new.env(parent = emptyenv()),
    now = now, max_attempts = max_attempts, cooldown = cooldown,
    global_max_attempts = global_max_attempts)
  state$attempts$wrong <- numeric()
  session$userData$rc_experimental <- state
  invisible(state)
}

#' Try an unlock code for the session.
#'
#' @return one of "unlocked", "wrong", "wait", "unavailable". Never the code.
#' @noRd
experimental_try_unlock <- function(session, code) {
  state <- session$userData$rc_experimental
  if (is.null(state) || !isTRUE(state$configured)) return("unavailable")
  t <- state$now()
  wrong <- state$attempts$wrong
  wrong <- wrong[wrong > t - 60]
  global <- .rc_experimental_env$wrong
  global <- if (is.null(global)) numeric() else global[global > t - 60]
  .rc_experimental_env$wrong <- global
  state$attempts$wrong <- wrong
  if (length(wrong) >= state$max_attempts ||
      (length(wrong) && t - max(wrong) < state$cooldown) ||
      length(global) >= state$global_max_attempts) {
    return("wait")
  }
  if (state$verify(code)) {
    state$attempts$wrong <- numeric()
    state$unlocked(rc_experimental_features())
    return("unlocked")
  }
  state$attempts$wrong <- c(wrong, t)
  .rc_experimental_env$wrong <- c(global, t)
  "wrong"
}

#' Is an experimental feature unlocked in this session?
#'
#' Generic gate for private features: FALSE outside a Shiny session, when no unlock
#' code is configured, or before the session typed the right code.
#' @param session Shiny session (module sessions share the root userData).
#' @param feature character, e.g. "explain".
#' @noRd
experimental_feature_enabled <- function(session = shiny::getDefaultReactiveDomain(), feature) {
  if (is.null(session)) return(FALSE)
  state <- tryCatch(session$userData$rc_experimental, error = function(e) NULL)
  if (is.null(state) || !isTRUE(state$configured)) return(FALSE)
  feature %in% state$unlocked()
}

#' Cache-key part for the session's experimental state ("" when locked).
#' @noRd
experimental_features_key <- function(session = shiny::getDefaultReactiveDomain()) {
  feats <- rc_experimental_features()
  on <- vapply(feats, function(f) experimental_feature_enabled(session, f), logical(1))
  paste(feats[on], collapse = ",")
}

#' Unlock box for the browser / observatory settings (NULL when not configured).
#' @noRd
experimental_unlock_ui <- function(ns, configured = experimental_unlock_configured()) {
  if (!isTRUE(configured)) return(NULL)
  fluidRow(
    column(4, passwordInput(ns("experimental_unlock_code"), "Experimental features (code)",
                            placeholder = "Unlock code")),
    column(2, tags$div(style = "margin-top: 25px;",
                       actionButton(ns("experimental_unlock_go"), "Unlock"))),
    column(6, tags$div(style = "margin-top: 30px;", textOutput(ns("experimental_unlock_status"),
                                                               inline = TRUE)))
  )
}

#' Settings rows of the experimental features (unlock box + server-rendered,
#' gated feature settings).
#' @noRd
experimental_settings_ui <- function(ns) {
  tagList(experimental_unlock_ui(ns), uiOutput(ns("explain_settings")))
}

#' Whether the app being built has an unlock configured (set by RiboCrypt_app).
#' @noRd
experimental_unlock_configured <- function() isTRUE(.rc_experimental_env$configured)

#' Module server part: unlock box handling and gated settings rendering.
#' @noRd
experimental_unlock_server <- function(input, output, session) {
  status <- shiny::reactiveVal("")
  output$experimental_unlock_status <- shiny::renderText(status())
  shiny::observeEvent(input$experimental_unlock_go, {
    code <- input$experimental_unlock_code
    res <- experimental_try_unlock(session, code)
    rm(code)
    shiny::updateTextInput(session, "experimental_unlock_code", value = "")
    status(switch(res, unlocked = "Unlocked", wait = "Too many attempts, try again in a minute",
                  "Wrong code"))
  }, ignoreInit = TRUE)
  output$explain_settings <- shiny::renderUI({
    shiny::req(experimental_feature_enabled(session, "explain"))
    explain_settings_ui(session$ns)
  })
  shiny::outputOptions(output, "explain_settings", suspendWhenHidden = FALSE)
  invisible(status)
}
