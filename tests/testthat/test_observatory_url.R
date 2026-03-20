test_that("observatory selection compaction and expansion roundtrip", {
  selections <- list(
    index = c("1", "2"),
    plot_selections = list(
      "1" = c("SRR1", "SRR2", "SRR3"),
      "2" = c("SRR2", "SRR4")
    ),
    data_table_selections = list(
      "1" = c("SRR1", "SRR2"),
      "2" = c("SRR4")
    ),
    labels = list("1" = "ctrl", "2" = "treated"),
    active_selection_id = "2"
  )

  compact <- RiboCrypt:::observatory_compact_selections(selections)
  expect_true(is.list(compact))
  expect_equal(compact$active, "2")
  expect_equal(compact$order, c("1", "2"))
  expect_true(length(compact$runs) >= 4)

  expanded <- RiboCrypt:::observatory_expand_selections(compact)
  expect_equal(expanded$index, c("1", "2"))
  expect_equal(expanded$active_selection_id, "2")
  expect_equal(expanded$labels[["1"]], "ctrl")
  expect_equal(expanded$labels[["2"]], "treated")
  expect_equal(sort(expanded$plot_selections[["1"]]), sort(c("SRR1", "SRR2", "SRR3")))
  expect_equal(sort(expanded$data_table_selections[["2"]]), "SRR4")
})

test_that("observatory state roundtrip via obs_state hash parameter", {
  state <- RiboCrypt:::observatory_state_from_inputs(
    selected_experiment = "all_samples-Homo_sapiens",
    color_by = c("tissue", "cell_line"),
    view = "browser",
    browser = list(
      gene = "ATF4-ENSG00000128272",
      tx = "ENST00000361436",
      frames_type = "area",
      kmer = 9,
      extendLeaders = 30,
      extendTrailers = 45,
      go = TRUE
    ),
    selections = list(
      index = c("1"),
      plot_selections = list("1" = c("SRR1", "SRR2")),
      data_table_selections = list("1" = c("SRR1")),
      labels = list("1" = "brain"),
      active_selection_id = "1"
    )
  )

  encoded <- RiboCrypt:::make_observatory_url_state_param(state)
  parsed <- RiboCrypt:::parse_observatory_url_hash(paste0("#Observatory?obs_state=", encoded))

  expect_equal(parsed$exp, "all_samples-Homo_sapiens")
  expect_equal(parsed$color_by, c("tissue", "cell_line"))
  expect_equal(parsed$view, "browser")
  expect_equal(parsed$browser$gene, "ATF4-ENSG00000128272")
  expect_equal(parsed$browser$tx, "ENST00000361436")
  expect_equal(parsed$browser$kmer, 9)
  expect_equal(parsed$selections$active_selection_id, "1")
  expect_equal(parsed$selections$labels[["1"]], "brain")
  expect_equal(sort(parsed$selections$plot_selections[["1"]]), c("SRR1", "SRR2"))
})

test_that("observatory parser supports legacy base64 json payload in query", {
  state <- list(exp = "legacy-exp", color_by = c("study"), view = "umap")
  legacy_encoded <- jsonlite::base64_enc(jsonlite::toJSON(state, auto_unbox = TRUE, null = "null"))
  parsed <- RiboCrypt:::parse_observatory_url_query(list(obs_state = legacy_encoded))

  expect_equal(parsed$exp, "legacy-exp")
  expect_equal(parsed$color_by, "study")
  expect_equal(parsed$view, "umap")
})

test_that("observatory query parsing is resilient to invalid payload", {
  expect_null(RiboCrypt:::parse_observatory_url_query(NULL))
  expect_null(RiboCrypt:::parse_observatory_url_query(list()))
  expect_null(RiboCrypt:::parse_observatory_url_query(list(obs_state = "not-valid-base64")))
  expect_null(RiboCrypt:::parse_observatory_url_hash("#Observatory"))
  expect_null(RiboCrypt:::parse_observatory_url_hash("#Observatory?obs_state=not-valid-base64"))
})

test_that("URL helpers can be explicitly targeted for observatory", {
  expect_no_error(RiboCrypt:::check_url_for_basic_parameters(mode = "observatory"))
  expect_no_error(RiboCrypt:::browser_specific_url_checker(target = "observatory"))
})

test_that("clipboard_url_text supports browser libraries and observatory run selections", {
  input <- list(
    dff = "exp-a",
    gene = "GENE1",
    tx = "TX1",
    library = c("lib-a", "lib-b"),
    unique_align = FALSE,
    frames_type = "columns",
    kmer = 1,
    log_scale = FALSE,
    log_scale_protein = FALSE,
    extendLeaders = 0,
    extendTrailers = 0,
    viewMode = FALSE,
    other_tx = FALSE,
    add_uorfs = FALSE,
    add_translon = FALSE,
    add_translons_transcode = FALSE,
    genomic_region = "",
    zoom_range = "",
    customSequence = "",
    phyloP = FALSE,
    mapability = FALSE,
    colors = "R",
    summary_track = FALSE,
    summary_track_type = "columns",
    collapsed_introns_width = 100,
    collapsed_introns = FALSE
  )
  session <- list(
    clientData = list(
      url_hostname = "localhost",
      url_port = "1234",
      url_pathname = "/app",
      url_hash = "#browser"
    )
  )

  browser_url <- RiboCrypt:::clipboard_url_text(
    input = input,
    session = session,
    mode = "browser",
    libraries = c("SRR100", "SRR200")
  )
  expect_match(browser_url, "library=SRR100,SRR200", fixed = TRUE)
  expect_match(browser_url, "#browser", fixed = TRUE)

  expected_state <- RiboCrypt:::observatory_state_from_inputs(
    selected_experiment = "obs-exp",
    color_by = c("study"),
    view = "browser",
    browser = list(
      gene = "GENE1",
      tx = "TX1",
      frames_type = "columns",
      kmer = 1,
      extendLeaders = 0,
      extendTrailers = 0,
      viewMode = FALSE,
      other_tx = FALSE,
      collapsed_introns = FALSE,
      collapsed_introns_width = 100,
      genomic_region = "",
      zoom_range = "",
      customSequence = "",
      go = TRUE
    ),
    selections = list(
      index = c("1"),
      plot_selections = list("1" = c("SRR100", "SRR200")),
      data_table_selections = list("1" = c("SRR100", "SRR200")),
      labels = list("1" = "brain"),
      active_selection_id = "1"
    )
  )

  observatory_url <- RiboCrypt:::clipboard_url_text(
    input = input,
    session = session,
    mode = "observatory",
    observatory = list(
      selected_experiment = function() "obs-exp",
      color_by = function() c("study"),
      selection_index = function() c("1"),
      library_selections = function() list("1" = c("SRR100", "SRR200")),
      library_selection_labels = function() list("1" = "brain"),
      active_selection_id = function() "1"
    )
  )
  expect_match(observatory_url, "/#Observatory\\?obs_state=")
  expect_equal(observatory_url, RiboCrypt:::make_observatory_url(expected_state, session))
})

test_that("observatory selection cache key is stable and sensitive to group membership", {
  selections <- list(
    "2" = c("SRR4", "SRR2", "SRR2"),
    "1" = c("SRR3", "SRR1")
  )
  labels <- list(
    "1" = "ctrl",
    "2" = "treated"
  )

  key <- RiboCrypt:::observatory_selection_cache_key(selections, labels)
  expect_equal(
    key,
    "2:treated:SRR2,SRR4|group|1:ctrl:SRR1,SRR3"
  )

  expect_equal(
    key,
    RiboCrypt:::observatory_selection_cache_key(
      list(
        "2" = c("SRR2", "SRR4"),
        "1" = c("SRR1", "SRR3")
      ),
      labels
    )
  )

  expect_false(identical(
    key,
    RiboCrypt:::observatory_selection_cache_key(
      list(
        "2" = c("SRR2", "SRR5"),
        "1" = c("SRR1", "SRR3")
      ),
      labels
    )
  ))
})
