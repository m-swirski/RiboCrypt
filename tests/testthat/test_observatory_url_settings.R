test_that("URL autostart waits for settings and the exact saved cohorts", {
  state <- list(view = "browser", browser = list(
    gene = "GENE", tx = "TX", kmer = 9, log_scale = TRUE,
    frames_subset = c("red", "blue"), go = TRUE),
    selections = list(index = "2", data_table_selections = list("2" = "SRR2")))
  input <- state$browser
  expect_true(observatory_browser_ready_to_kickoff(state, input, list("2" = "SRR2")))
  expect_false(observatory_browser_ready_to_kickoff(state, input, list("1" = "SRR2")))
  expect_false(observatory_browser_ready_to_kickoff(state, input, list("2" = c("SRR1", "SRR2"))))
  for (field in c("kmer", "log_scale", "frames_subset")) {
    pending <- input
    pending[[field]] <- NULL
    expect_false(observatory_browser_ready_to_kickoff(state, pending, list("2" = "SRR2")))
  }
  state$browser$go <- FALSE
  expect_false(observatory_browser_ready_to_kickoff(state, input, list("2" = "SRR2")))
})

test_that("URL control restoration sends all supported values including false and empty", {
  values <- list(gene = "GENE", tx = "TX", kmer = 9, log_scale = TRUE,
    log_scale_protein = FALSE, phyloP = FALSE, mapability = TRUE,
    colors = "Color_blind", summary_track = TRUE, summary_track_type = "area",
    withFrames = FALSE, frames_subset = character(), customSequence = "",
    extendLeaders = 10, collapsed_introns = TRUE)
  messages <- list()
  shiny::testServer(function(input, output, session) {
    session$sendInputMessage <- function(inputId, message) messages[[inputId]] <<- message
    apply_observatory_browser_url_state(session, values)
  }, {
    expect_setequal(names(messages), names(values))
    for (field in names(values)) {
      expect_equal(as.character(messages[[field]]$value), as.character(values[[field]]), info = field)
    }
  })
})

test_that("malformed Observatory payloads are rejected without breaking startup", {
  valid <- observatory_state_from_inputs("exp", "tissue", "browser",
    list(gene = "G", tx = "T", go = TRUE),
    list(index = "2", active_selection_id = "2", labels = list("2" = "Saved"),
      plot_selections = list("2" = "SRR2"), data_table_selections = list("2" = "SRR2")))
  expect_true(is.list(parse_observatory_url_state_param(make_observatory_url_state_param(valid))))
  malformed <- list(list(v = 2), list(browser = list(kmer = c(1, 2))),
    list(browser = list(log_scale = "invalid")), list(browser = list(colors = "invalid")),
    list(selections = list(order = c("1", "1"), active = "1")),
    list(selections = list(order = "1", active = "1", runs = "SRR1", d = list("1" = 2))))
  for (state in malformed) {
    expect_null(parse_observatory_url_state_param(make_observatory_url_state_param(state)))
  }
  expect_null(parse_observatory_url_state_param(NA_character_))
  expect_null(parse_observatory_url_state_param(c("a", "b")))
})

test_that("annotation fallbacks also update effective URL settings", {
  names <- data.table::data.table(value = c("T1", "T2"), label = c("G", "G"))
  options <- c(default_gene_meta = "G", default_isoform_meta = "T1")
  for (tx in c("", "missing", "T2")) {
    state <- list(view = "browser", browser = list(gene = "G", tx = tx, go = TRUE))
    defaults <- observatory_browser_url_defaults(options, state, names)
    resolved <- observatory_resolve_browser_url_state(state, defaults)
    expect_equal(resolved$browser$tx, if (tx == "T2") "T2" else "T1")
    expect_true(observatory_browser_ready_to_kickoff(resolved, resolved$browser, list("1" = "SRR1")))
  }
})

test_that("shared URLs preserve HTTPS and reverse proxy paths", {
  session <- list(clientData = list(url_protocol = "https:", url_hostname = "example.org",
    url_port = "", url_pathname = "/apps/ribo/"))
  expect_equal(getHostFromURL(session), "https://example.org/apps/ribo")
  session$clientData$url_port <- "8443"
  expect_equal(getHostFromURL(session), "https://example.org:8443/apps/ribo")
  expect_equal(getPageFromURL(url = "#Observatory?obs_state=abc"), "Observatory")
})

test_that("URL collection annotations are loaded together and only when needed", {
  calls <- character()
  local_mocked_bindings(
    name = function(df) df$id,
    get_exp = function(experiment, ...) list(id = experiment),
    get_gene_name_categories = function(df) paste0(df$id, "-names"),
    loadRegion = function(df, region = "mrna") {
      calls <<- c(calls, region)
      paste(df$id, region)
    }, .package = "RiboCrypt")
  initial <- list(id = "human")
  for (state in list(NULL, list(exp = "human"), list(exp = "unknown"))) {
    expect_null(observatory_url_collection_init(state, c("human", "yeast"), initial, ""))
  }
  expect_length(calls, 0)
  restored <- observatory_url_collection_init(list(exp = "yeast"), c("human", "yeast"), initial, "")
  expect_equal(restored, list(df = list(id = "yeast"), names = "yeast-names",
    tx = "yeast mrna", cds = "yeast cds"))
  expect_equal(calls, c("mrna", "cds"))
})
