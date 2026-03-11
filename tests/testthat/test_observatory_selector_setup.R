make_observatory_selector_fixture <- function() {
  all_exp <- data.frame(
    name = c("exp-a", "exp-b"),
    stringsAsFactors = FALSE
  )

  libraries_df <- data.table::data.table(
    Run = c("SRR1", "SRR2", "SRR3", "SRR4", "SRR5"),
    LIBRARYTYPE = c("RFP", "RFP", "RFP", "RFP", "RFP"),
    TISSUE = c("Cervix", "Cervix", "Lung", "Cervix", "Cervix"),
    CELL_LINE = c("HeLa", "HeLa", "A549", "SiHa", "HeLa"),
    BioProject = c("PRJ1", "PRJ1", "PRJ2", "PRJ3", "PRJ4"),
    author = c("Auth1", "Auth1", "Auth2", "Auth3", "Auth4"),
    inhibitors = c("none", "none", "drug", "none", "drug"),
    tissues_cell_lines = c(
      "Cervix | HeLa",
      "Cervix | HeLa",
      "Lung | A549",
      "Cervix | SiHa",
      "Cervix | HeLa"
    )
  )

  umap_df <- data.table::copy(libraries_df)[, `:=`(
    sample = Run,
    color_column = TISSUE,
    `UMAP 1` = c(0, 1, 2, 3, 4),
    `UMAP 2` = c(4, 3, 2, 1, 0)
  )]
  data.table::setattr(umap_df, "color.by", "tissue")

  experiment_df <- ORFik::ORFik.template.experiment()[9:10, ]

  list(
    all_exp = all_exp,
    libraries_df = libraries_df,
    umap_df = umap_df,
    experiment_df = experiment_df,
    browser_options = c(default_experiment_meta = "exp-a")
  )
}

observatory_selector_harness_server <- function(
  id,
  all_exp,
  experiment_df,
  libraries_df,
  browser_options = c(),
  initial_url_state = NULL
) {
  shiny::moduleServer(id, function(input, output, session) {
    observatory_url_state <- shiny::reactiveVal(initial_url_state)
    selections <- RiboCrypt:::observatory_selector_server(
      "selector",
      all_exp,
      shiny::reactive(experiment_df),
      libraries_df,
      experiments = NULL,
      org = shiny::reactive(NULL),
      rv = shiny::reactiveValues(),
      browser_options = browser_options,
      observatory_url_state = observatory_url_state
    )

    selected_libraries <- selections$selected_libraries
    selected_experiment <- selections$selected_experiment
    color_by <- selections$color_by
    active_selection_id <- selections$active_selection_id
  })
}

test_that("observatory_ui includes selector and browser tabs", {
  ui <- RiboCrypt:::observatory_ui(
    "observatory",
    data.frame(name = "exp-a", stringsAsFactors = FALSE),
    c(
      default_experiment_meta = "exp-a",
      default_frame_type = "area",
      default_kmer = "1"
    )
  )

  html <- as.character(ui)
  expect_match(html, "Observatory")
  expect_match(html, "Select libraries")
  expect_match(html, "Browse")
})

test_that("observatory selector harness initializes URL-backed selection state", {
  fixture <- make_observatory_selector_fixture()
  initial_state <- list(
    exp = "exp-a",
    color_by = c("tissue", "cell_line"),
    selections = list(
      index = c("1", "2"),
      plot_selections = list(
        "1" = c("SRR1", "SRR2"),
        "2" = c("SRR3")
      ),
      data_table_selections = list(
        "1" = c("SRR1"),
        "2" = c("SRR3")
      ),
      labels = list(
        "1" = "brain",
        "2" = "heart"
      ),
      active_selection_id = "2"
    )
  )

  local_mocked_bindings(
    allsamples_observer_controller = function(input, output, session) invisible(NULL),
    create_observatory_module = function(meta_experiment_df, libraries_df) {
      list(
        get_libraries_data = function(library_types = c("RFP")) {
          libraries_df[LIBRARYTYPE %in% library_types]
        },
        get_umap_data = function(color_by = c("tissue", "cell_line")) {
          data.table::copy(fixture$umap_df)
        }
      )
    },
    .package = "RiboCrypt"
  )

  shiny::testServer(
    observatory_selector_harness_server,
    args = list(
      all_exp = fixture$all_exp,
      experiment_df = fixture$experiment_df,
      libraries_df = fixture$libraries_df,
      browser_options = fixture$browser_options,
      initial_url_state = initial_state
    ),
    {
      session$setInputs(
        `selector-dff` = "exp-a",
        `selector-color_by` = c("tissue", "cell_line")
      )
      session$flushReact()

      expect_equal(selected_experiment(), "exp-a")
      expect_equal(color_by(), c("tissue", "cell_line"))
      expect_equal(active_selection_id(), "2")
      expect_equal(selected_libraries$index(), c("1", "2"))
      expect_equal(selected_libraries$labels()[["1"]], "brain")
    }
  )
})

test_that("observatory selector derives data-table selection from UMAP selection", {
  fixture <- make_observatory_selector_fixture()

  local_mocked_bindings(
    allsamples_observer_controller = function(input, output, session) invisible(NULL),
    create_observatory_module = function(meta_experiment_df, libraries_df) {
      list(
        get_libraries_data = function(library_types = c("RFP")) {
          libraries_df[LIBRARYTYPE %in% library_types]
        },
        get_umap_data = function(color_by = c("tissue", "cell_line")) {
          data.table::copy(fixture$umap_df)
        }
      )
    },
    .package = "RiboCrypt"
  )

  shiny::testServer(
    observatory_selector_harness_server,
    args = list(
      all_exp = fixture$all_exp,
      experiment_df = fixture$experiment_df,
      libraries_df = fixture$libraries_df,
      browser_options = fixture$browser_options
    ),
    {
      session$setInputs(
        `selector-dff` = "exp-a",
        `selector-color_by` = c("tissue"),
        `selector-go` = 1
      )
      session$flushReact()

      session$setInputs(`selector-libraries_umap_plot_selection` = c("SRR1", "SRR3"))
      session$flushReact()

      expect_equal(active_selection_id(), "1")
      expect_equal(
        sort(selected_libraries$data_table_selections()[["1"]]),
        c("SRR1", "SRR3")
      )
    }
  )
})

test_that("observatory selector auto-labels UMAP selections by study first", {
  fixture <- make_observatory_selector_fixture()

  local_mocked_bindings(
    allsamples_observer_controller = function(input, output, session) invisible(NULL),
    create_observatory_module = function(meta_experiment_df, libraries_df) {
      list(
        get_libraries_data = function(library_types = c("RFP")) {
          libraries_df[LIBRARYTYPE %in% library_types]
        },
        get_umap_data = function(color_by = c("tissue", "cell_line")) {
          data.table::copy(fixture$umap_df)
        }
      )
    },
    .package = "RiboCrypt"
  )

  shiny::testServer(
    observatory_selector_harness_server,
    args = list(
      all_exp = fixture$all_exp,
      experiment_df = fixture$experiment_df,
      libraries_df = fixture$libraries_df,
      browser_options = fixture$browser_options
    ),
    {
      session$setInputs(`selector-dff` = "exp-a", `selector-color_by` = c("tissue"), `selector-go` = 1)
      session$flushReact()

      session$setInputs(`selector-libraries_umap_plot_selection` = c("SRR1", "SRR2"))
      session$flushReact()

      expect_equal(selected_libraries$labels()[["1"]], "PRJ1 subset")
    }
  )
})

test_that("observatory selector auto-labels UMAP selections by cell line when study is mixed", {
  fixture <- make_observatory_selector_fixture()

  local_mocked_bindings(
    allsamples_observer_controller = function(input, output, session) invisible(NULL),
    create_observatory_module = function(meta_experiment_df, libraries_df) {
      list(
        get_libraries_data = function(library_types = c("RFP")) {
          libraries_df[LIBRARYTYPE %in% library_types]
        },
        get_umap_data = function(color_by = c("tissue", "cell_line")) {
          data.table::copy(fixture$umap_df)
        }
      )
    },
    .package = "RiboCrypt"
  )

  shiny::testServer(
    observatory_selector_harness_server,
    args = list(
      all_exp = fixture$all_exp,
      experiment_df = fixture$experiment_df,
      libraries_df = fixture$libraries_df,
      browser_options = fixture$browser_options
    ),
    {
      session$setInputs(`selector-dff` = "exp-a", `selector-color_by` = c("tissue"), `selector-go` = 1)
      session$flushReact()

      session$setInputs(`selector-libraries_umap_plot_selection` = c("SRR1", "SRR5"))
      session$flushReact()

      expect_equal(selected_libraries$labels()[["1"]], "HeLa subset")
    }
  )
})

test_that("observatory selector auto-labels UMAP selections by tissue when study and cell line are mixed", {
  fixture <- make_observatory_selector_fixture()

  local_mocked_bindings(
    allsamples_observer_controller = function(input, output, session) invisible(NULL),
    create_observatory_module = function(meta_experiment_df, libraries_df) {
      list(
        get_libraries_data = function(library_types = c("RFP")) {
          libraries_df[LIBRARYTYPE %in% library_types]
        },
        get_umap_data = function(color_by = c("tissue", "cell_line")) {
          data.table::copy(fixture$umap_df)
        }
      )
    },
    .package = "RiboCrypt"
  )

  shiny::testServer(
    observatory_selector_harness_server,
    args = list(
      all_exp = fixture$all_exp,
      experiment_df = fixture$experiment_df,
      libraries_df = fixture$libraries_df,
      browser_options = fixture$browser_options
    ),
    {
      session$setInputs(`selector-dff` = "exp-a", `selector-color_by` = c("tissue"), `selector-go` = 1)
      session$flushReact()

      session$setInputs(`selector-libraries_umap_plot_selection` = c("SRR4", "SRR5"))
      session$flushReact()

      expect_equal(selected_libraries$labels()[["1"]], "Cervix subset")
    }
  )
})

test_that("observatory selector respects explicit DT-selected run ids", {
  fixture <- make_observatory_selector_fixture()

  local_mocked_bindings(
    allsamples_observer_controller = function(input, output, session) invisible(NULL),
    create_observatory_module = function(meta_experiment_df, libraries_df) {
      list(
        get_libraries_data = function(library_types = c("RFP")) {
          libraries_df[LIBRARYTYPE %in% library_types]
        },
        get_umap_data = function(color_by = c("tissue", "cell_line")) {
          data.table::copy(fixture$umap_df)
        }
      )
    },
    .package = "RiboCrypt"
  )

  shiny::testServer(
    observatory_selector_harness_server,
    args = list(
      all_exp = fixture$all_exp,
      experiment_df = fixture$experiment_df,
      libraries_df = fixture$libraries_df,
      browser_options = fixture$browser_options
    ),
    {
      session$setInputs(
        `selector-dff` = "exp-a",
        `selector-color_by` = c("tissue"),
        `selector-go` = 1
      )
      session$flushReact()

      session$setInputs(`selector-libraries_data_table_selected_runs` = c("SRR2"))
      session$flushReact()

      expect_equal(
        selected_libraries$data_table_selections()[["1"]],
        "SRR2"
      )
    }
  )
})

test_that("observatory selector derives active selection label from DT search terms", {
  fixture <- make_observatory_selector_fixture()

  local_mocked_bindings(
    allsamples_observer_controller = function(input, output, session) invisible(NULL),
    create_observatory_module = function(meta_experiment_df, libraries_df) {
      list(
        get_libraries_data = function(library_types = c("RFP")) {
          libraries_df[LIBRARYTYPE %in% library_types]
        },
        get_umap_data = function(color_by = c("tissue", "cell_line")) {
          data.table::copy(fixture$umap_df)
        }
      )
    },
    .package = "RiboCrypt"
  )

  shiny::testServer(
    observatory_selector_harness_server,
    args = list(
      all_exp = fixture$all_exp,
      experiment_df = fixture$experiment_df,
      libraries_df = fixture$libraries_df,
      browser_options = fixture$browser_options
    ),
    {
      session$setInputs(
        `selector-dff` = "exp-a",
        `selector-color_by` = c("tissue"),
        `selector-go` = 1
      )
      session$flushReact()

      session$setInputs(`selector-libraries_data_table_manual_search` = "Brain")
      session$flushReact()

      expect_match(selected_libraries$labels()[["1"]], "brain")
    }
  )
})
