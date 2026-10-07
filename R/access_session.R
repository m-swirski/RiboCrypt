#' Configure Authenticated RiboCrypt Sessions
#' @param database Existing access-registry SQLite file (not `:memory:`).
#' @param issuer Fixed OIDC issuer URL.
#' @param gateway_secret Shared secret injected by the authentication gateway.
#' @param login_url Same-site login route.
#' @param logout_url Same-site route that clears proxy and identity-provider sessions.
#' @param refresh_ms Permission-revocation polling interval in milliseconds.
#' @return Configuration for the `access_control` argument of `RiboCrypt_app()`.
#' @export
ribocrypt_access_control <- function(database, issuer, gateway_secret,
                                    login_url = "/oauth2/start?rd=/",
                                    logout_url = "/ribocrypt/logout", refresh_ms = 1000L) {
  if (!file.exists(database)) stop("Create the access registry before starting the app.", call. = FALSE)
  ribocrypt_access_identity(list(HTTP_RIBOCRYPT_GATEWAY_SECRET = gateway_secret), issuer, gateway_secret)
  for (url in c(login_url, logout_url))
    if (!grepl("^/[^/]", url) || grepl("[[:cntrl:]\\\\]", url)) stop("Account routes must be same-site paths.")
  if (length(refresh_ms) != 1L || !is.numeric(refresh_ms) || !is.finite(refresh_ms) || refresh_ms < 100)
    stop("refresh_ms must be a finite number of at least 100.", call. = FALSE)
  structure(list(database = normalizePath(database), issuer = issuer, gateway_secret = gateway_secret,
                 login_url = login_url, logout_url = logout_url, refresh_ms = refresh_ms),
            class = "ribocrypt_access_control")
}

#' @rdname ribocrypt_access_control
#' @param x Access-control configuration to print, with its secret redacted.
#' @param ... Additional arguments (unused).
#' @export
print.ribocrypt_access_control <- function(x, ...) {
  cat("RiboCrypt access control\nIssuer:", x$issuer,
      "\nDatabase:", x$database, "\nGateway secret: [redacted]\n")
  invisible(x)
}

rc_access_context <- function() {
  session <- shiny::getDefaultReactiveDomain()
  if (is.null(session)) return(NULL)
  session$userData$ribocrypt_access
}

rc_access_check <- function(context = rc_access_context(), capability = "read") {
  if (is.null(context)) return(invisible(TRUE))
  current <- ribocrypt_access_catalog(context$con, context$identity, capability)
  if (!all(context$catalog$id %in% current$id)) stop("Dataset permissions changed or access denied. Reload the page.", call. = FALSE)
  invisible(TRUE)
}

rc_access_export_buttons <- function(buttons) {
  allowed <- tryCatch({ rc_access_check(capability = "download"); TRUE }, error = function(e) FALSE)
  if (allowed) buttons else list()
}

rc_access_account_tab <- function() {
  context <- rc_access_context()
  if (is.null(context)) return(NULL)
  shiny::tabPanel("Account", value = "Account", icon = shiny::icon("user"),
    shiny::h3("Account"), shiny::p(if (is.null(context$identity)) "Anonymous" else context$identity$subject),
    shiny::h4("Assigned datasets"), DT::datatable(context$catalog[, c("id", "experiment", "public")],
      rownames = FALSE, options = list(pageLength = 15, dom = "t")))
}

rc_access_heartbeat_script <- function(session) {
  url <- session$registerDataObj("rc_account_check", NULL, function(data, req) {
    list(status = 204L, headers = list(`Cache-Control` = "no-store"), body = "")
  })
  shiny::tags$script(shiny::HTML(sprintf(paste0(
    "(function() { const check = function() { fetch(%s, {credentials:'same-origin',cache:'no-store'})",
    ".then(function(response) { if (!response.ok) { document.getElementById('authorized_app').replaceChildren();",
    "Shiny.shinyapp.$socket.close(); window.location.reload(); } }).catch(function() {}); };",
    "setInterval(check, 30000); document.addEventListener('visibilitychange', function() { if (!document.hidden) check(); }); })();"),
    jsonlite::toJSON(url, auto_unbox = TRUE))))
}

rc_read_experiment <- function(file, in.dir = ORFik::config()["exp"],
                               validate = FALSE, output.env = .GlobalEnv) {
  context <- rc_access_context()
  if (!is.null(context)) {
    rc_access_check(context)
    row <- context$catalog[context$catalog$experiment == file, , drop = FALSE]
    if (nrow(row) != 1L) stop("Dataset unavailable or access denied.", call. = FALSE)
    in.dir <- row$directory
    file <- normalizePath(file.path(in.dir, paste0(row$experiment, ".csv")), mustWork = TRUE)
    if (identical(output.env, .GlobalEnv)) output.env <- context$env
  }
  df <- ORFik::read.experiment(file, in.dir = in.dir, validate = validate, output.env = output.env)
  if (!is.null(context) && !identical(ORFik::name(df), row$experiment))
    stop("Registered experiment name must match the CSV name field.", call. = FALSE)
  df
}

rc_access_metadata <- function(metadata, context) {
  if (is.null(metadata)) return(NULL)
  if (is.character(metadata)) metadata <- data.table::fread(metadata)
  if (!"Run" %in% names(metadata)) stop("Metadata must contain Run.", call. = FALSE)
  data.table::copy(data.table::as.data.table(metadata)[Run %in% context$runs])
}

rc_access_cache <- function(cache, context) {
  guarded <- cache
  guarded$get <- function(key) { rc_access_check(context); cache$get(key) }
  guarded$set <- function(key, value) { rc_access_check(context); cache$set(key, value) }
  guarded
}

rc_download_handler <- function(filename, content, ...) {
  downloadHandler(filename = filename, content = function(file) {
    rc_access_check(capability = "download")
    content(file)
  }, ...)
}

rc_access_response <- function(status = 403L) {
  list(status = status, headers = list(`Content-Type` = "text/plain", `Cache-Control` = "no-store"),
       body = "Access denied. Reload the page or sign in.")
}

rc_access_request_handler <- function(handler, context, config) {
  force(handler)
  function(req) {
    allowed <- tryCatch({
      identity <- ribocrypt_access_identity(req, config$issuer, config$gateway_secret)
      if (!identical(identity, context$identity)) stop("Session identity mismatch.")
      rc_access_check(context, if (grepl("^/(download|file)/", req$PATH_INFO)) "download" else "read")
      TRUE
    }, error = function(e) FALSE)
    if (!allowed) return(rc_access_response())
    handler(req)
  }
}

rc_access_attach <- function(session, config, request = session$request) {
  identity <- ribocrypt_access_identity(request, config$issuer, config$gateway_secret)
  con <- ribocrypt_access_db(config$database)
  session$onSessionEnded(function() DBI::dbDisconnect(con))
  context <- list(con = con, identity = identity, catalog = ribocrypt_access_catalog(con, identity), env = new.env())
  if (anyDuplicated(context$catalog$experiment)) stop("Authorized experiments need distinct names in this deployment.")
  rc_access_record_account(con, identity)
  session$userData$ribocrypt_access <- context
  shiny::shinyOptions(cache = rc_access_cache(session$cache, context))
  rc_access_protect_session_http(session, context, config)
  shiny::observe({
    shiny::invalidateLater(config$refresh_ms, session)
    tryCatch(rc_access_check(context), error = function(e) {
      session$cache$reset()
      session$sendCustomMessage("ribocrypt-access-revoked", list())
      session$close()
    })
  }, priority = 1000)
  context
}

# Shiny routes session downloads/data objects before the app HTTP handler.
rc_access_protect_session_http <- function(session, context, config) {
  if (!is.function(session$handleRequest)) return(invisible(NULL))
  handler <- rc_access_request_handler(session$handleRequest, context, config)
  locked <- bindingIsLocked("handleRequest", session)
  if (locked) unlockBinding("handleRequest", session)
  on.exit(if (locked) lockBinding("handleRequest", session))
  session$handleRequest <- handler
  invisible(NULL)
}

rc_access_record_account <- function(con, identity) {
  if (is.null(identity)) return(invisible(NULL))
  DBI::dbExecute(con, paste0("INSERT INTO accounts (issuer,subject,first_seen,last_seen) VALUES (?,?,?,?) ",
    "ON CONFLICT(issuer,subject) DO UPDATE SET last_seen=excluded.last_seen"),
    params = list(identity$issuer, identity$subject, as.numeric(Sys.time()), as.numeric(Sys.time())))
  invisible(identity)
}

rc_access_experiment_catalog <- function(context) {
  experiments <- lapply(seq_len(nrow(context$catalog)), function(i) {
    row <- context$catalog[i, ]
    df <- rc_read_experiment(row$experiment)
    list(df = df, info = data.table::data.table(name = row$experiment, organism = ORFik::organism(df),
      author = df@author, libtypes = list(unique(df$libtype)), samples = nrow(df)))
  })
  context$runs <- unique(unlist(lapply(experiments, function(x) ORFik::runIDs(x$df))))
  context$runs <- context$runs[!is.na(context$runs) & nzchar(context$runs)]
  context$collection_roots <- unique(unlist(lapply(experiments, function(x)
    c(collection_dir_from_exp(x$df, new_format = TRUE), collection_dir_from_exp(x$df, new_format = FALSE)))))
  list(context = context, all_exp = data.table::rbindlist(lapply(experiments, `[[`, "info")))
}

rc_access_defaults <- function(options, experiments, collections) {
  for (suffix in c("", "_meta")) {
    key <- paste0("default_experiment", suffix)
    allowed <- if (suffix == "") experiments else collections
    if (!length(options[key]) || is.na(options[key]) || !options[key] %in% allowed) {
      options <- options[!names(options) %in% paste0(c("default_experiment", "default_gene", "default_isoform"), suffix)]
      if (length(allowed)) options[key] <- allowed[1]
      if (suffix == "") options <- options[!names(options) %in% "default_libs"]
    }
  }
  if (!is.na(options["default_experiment_translon"]) && !options["default_experiment_translon"] %in% experiments)
    options <- options[!names(options) %in% "default_experiment_translon"]
  options
}

rc_access_account_ui <- function(context, config) {
  logged_in <- !is.null(context$identity)
  shiny::div(style = "display:flex;gap:16px;align-items:center;padding:8px 16px;border-bottom:1px solid #ddd;",
    shiny::span(if (logged_in) "Signed in" else "Public access"),
    shiny::tags$a(href = if (logged_in) config$logout_url else config$login_url,
      if (logged_in) "Sign out" else "Log in / Create account"))
}

rc_access_ready_script <- function() {
  shiny::tags$script(shiny::HTML(paste0(
    "$(document).on('shiny:bound.rcAccess', function(event) {",
    "if (event.target.id !== 'navbarID') return;",
    "$(document).off('shiny:bound.rcAccess');",
    "setTimeout(function() { Shiny.setInputValue('rc_authorized_ready', true, {priority:'event'}); }, 0);",
    "});")))
}

rc_access_app <- function(config, args, build = RiboCrypt_app) {
  stopifnot(inherits(config, "ribocrypt_access_control"))
  ui <- bslib::page_fluid(theme = rc_theme(), shinyjs::useShinyjs(), rclipboard::rclipboardSetup(),
    shiny::uiOutput("authorized_app"),
    shiny::tags$script(shiny::HTML("Shiny.addCustomMessageHandler('ribocrypt-access-revoked', function(message) { document.getElementById('authorized_app').replaceChildren(); window.location.reload(); });")))
  server <- function(input, output, session) {
    context <- rc_access_attach(session, config)
    if (!nrow(context$catalog)) {
      output$authorized_app <- shiny::renderUI(shiny::tagList(rc_access_account_ui(context, config),
        shiny::p("No datasets assigned to this account.")))
      return()
    }
    prepared <- rc_access_experiment_catalog(context)
    context <- prepared$context
    session$userData$ribocrypt_access <- context
    args$all_exp <- prepared$all_exp
    collection_names <- if (is.null(args$all_exp_meta)) grep("all_samples-", args$all_exp$name, value = TRUE) else args$all_exp_meta$name
    args$all_exp_meta <- args$all_exp[name %in% collection_names]
    args$metadata <- rc_access_metadata(args$metadata, context)
    args$browser_options <- rc_access_defaults(args$browser_options, args$all_exp$name, args$all_exp_meta$name)
    child <- do.call(build, args)
    heartbeat <- rc_access_heartbeat_script(session)
    output$authorized_app <- shiny::renderUI(shiny::tagList(rc_access_account_ui(context, config),
      attr(child, "ribocrypt_ui"), rc_access_ready_script(), heartbeat))
    shiny::observeEvent(input$rc_authorized_ready, {
      child$serverFuncSource()(input, output, session)
    }, once = TRUE, ignoreInit = FALSE, priority = 100)
  }
  app <- shiny::shinyApp(ui, server, options = args$options)
  handler <- app$httpHandler
  app$httpHandler <- function(req) {
    valid <- tryCatch({ ribocrypt_access_identity(req, config$issuer, config$gateway_secret); TRUE }, error = function(e) FALSE)
    if (!valid) return(rc_access_response())
    handler(req)
  }
  app
}
