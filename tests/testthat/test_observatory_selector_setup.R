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

test_that("observatory UMAP JS attaches robust double-click reset handlers", {
  render_code <- RiboCrypt:::fetchJS("umap_plot_extension.js")

  expect_match(render_code, "plotly_doubleclick")
  expect_match(render_code, "plotly_deselect")
  expect_match(render_code, "addEventListener\\(\"dblclick\"")
  expect_match(render_code, "querySelectorAll")
  expect_match(render_code, "__rcUmapDblclickBound")
  expect_match(render_code, "preventDefault")
  expect_match(render_code, "stopImmediatePropagation")
  expect_match(render_code, "sendSelection\\(\\[\\]\\)")
  expect_match(render_code, "selectedpoints\"\\] = null")
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

test_that("observatory selector defaults active selection label to All merged", {
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

      expect_equal(selected_libraries$labels()[["1"]], "All merged")
      expect_equal(sort(selected_libraries$data_table_selections()[["1"]]), sort(fixture$libraries_df$Run))
    }
  )
})

test_that("observatory selection helper treats full-run selection as All merged", {
  expect_true(RiboCrypt:::observatory_selection_is_all_merged(NULL, c("SRR1", "SRR2")))
  expect_true(RiboCrypt:::observatory_selection_is_all_merged(c("SRR2", "SRR1"), c("SRR1", "SRR2")))
  expect_false(RiboCrypt:::observatory_selection_is_all_merged(c("SRR1"), c("SRR1", "SRR2")))
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

test_that("observatory table prioritizes key metadata columns", {
  df <- data.table::data.table(
    score = c(1.2345, 6.789),
    author = c("Auth1", "Auth2"),
    Run = c("SRR1", "SRR2"),
    BioProject = c("PRJ1", "PRJ2"),
    extra = c("x", "y"),
    TISSUE = c("Liver", "Brain")
  )

  expect_equal(
    RiboCrypt:::observatory_priority_library_columns(df),
    c("Run", "BioProject", "TISSUE", "author", "score", "extra")
  )
})

test_that("observatory table hides library type and scientific name columns", {
  df <- data.table::data.table(
    Run = c("SRR1", "SRR2"),
    BioProject = c("PRJ1", "PRJ2"),
    LIBRARYTYPE = c("RFP", "RFP"),
    ScientificName = c("Homo sapiens", "Homo sapiens"),
    author = c("Auth1", "Auth2")
  )

  formatted <- RiboCrypt:::observatory_format_libraries_df(df)

  expect_false("LIBRARYTYPE" %in% colnames(formatted))
  expect_false("ScientificName" %in% colnames(formatted))
  expect_equal(colnames(formatted), c("Run", "BioProject", "author"))
})

test_that("observatory selector applies DT column filters using displayed column order", {
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

      displayed_columns <- colnames(RiboCrypt:::observatory_format_libraries_df(fixture$libraries_df))
      column_filters <- rep("", length(displayed_columns))
      column_filters[match("BioProject", displayed_columns)] <- "PRJ1"

      session$setInputs(`selector-libraries_data_table_manual_search_columns` = column_filters)
      session$flushReact()

      expect_equal(
        sort(selected_libraries$data_table_selections()[["1"]]),
        c("SRR1", "SRR2")
      )
    }
  )
})

test_that("observatory selector clears active data-table selection when plot selection resets", {
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
      expect_equal(sort(selected_libraries$data_table_selections()[["1"]]), c("SRR1", "SRR3"))

      session$setInputs(`selector-libraries_umap_plot_selection` = character())
      session$flushReact()
      expect_equal(
        sort(selected_libraries$data_table_selections()[["1"]]),
        sort(fixture$libraries_df$Run)
      )
      expect_equal(selected_libraries$labels()[["1"]], "All merged")
    }
  )
})

test_that("observatory selector clears stale plot subset when go is triggered again", {
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
      expect_equal(sort(selected_libraries$data_table_selections()[["1"]]), c("SRR1", "SRR3"))

      session$setInputs(`selector-go` = 2)
      session$flushReact()

      expect_equal(sort(selected_libraries$plot_selections()[["1"]]), sort(fixture$libraries_df$Run))
      expect_equal(sort(selected_libraries$data_table_selections()[["1"]]), sort(fixture$libraries_df$Run))
      expect_equal(selected_libraries$labels()[["1"]], "All merged")
    }
  )
})

test_that("observatory selector new selections start as All merged over all runs", {
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

      session$setInputs(`selector-library_selection-active_selection_id` = "New selection...")
      session$flushReact()

      expect_equal(active_selection_id(), "2")
      expect_equal(selected_libraries$labels()[["2"]], "All merged")
      expect_equal(sort(selected_libraries$plot_selections()[["2"]]), sort(fixture$libraries_df$Run))
      expect_equal(sort(selected_libraries$data_table_selections()[["2"]]), sort(fixture$libraries_df$Run))
    }
  )
})

test_that("observatory selector resets active subset label to All merged on plot reset", {
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

      session$setInputs(`selector-library_selection-active_selection_id` = "New selection...")
      session$flushReact()
      expect_equal(active_selection_id(), "2")

      session$setInputs(`selector-libraries_umap_plot_selection` = c("SRR1", "SRR3"))
      session$flushReact()

      session$setInputs(`selector-libraries_umap_plot_selection` = character())
      session$flushReact()
      expect_equal(selected_libraries$labels()[["2"]], "All merged")
    }
  )
})

test_that("observatory selector resets active subset label to All merged on NULL plot reset", {
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

      session$setInputs(`selector-library_selection-active_selection_id` = "New selection...")
      session$flushReact()
      expect_equal(active_selection_id(), "2")

      session$setInputs(`selector-libraries_umap_plot_selection` = c("SRR1", "SRR3"))
      session$flushReact()

      session$setInputs(`selector-libraries_umap_plot_selection` = NULL)
      session$flushReact()
      expect_equal(selected_libraries$labels()[["2"]], "All merged")
    }
  )
})

test_that("observatory browser resolves empty active selection to all runs", {
  selections <- list("1" = NULL, "2" = c("SRR9"))

  resolved <- RiboCrypt:::observatory_resolve_library_selections(
    library_selections = selections,
    all_library_runs = c("SRR1", "SRR2", "SRR3"),
    active_selection_id = "1",
    library_selection_labels = list("1" = "All merged", "2" = "subset")
  )

  expect_equal(resolved[["1"]], c("SRR1", "SRR2", "SRR3"))
  expect_equal(resolved[["2"]], "SRR9")
})

test_that("observatory browser creates default all-runs selection when none exists", {
  resolved <- RiboCrypt:::observatory_resolve_library_selections(
    library_selections = list(),
    all_library_runs = c("SRR1", "SRR2"),
    active_selection_id = "5",
    library_selection_labels = list()
  )

  expect_equal(names(resolved), "5")
  expect_equal(resolved[["5"]], c("SRR1", "SRR2"))
})

test_that("observatory browser preserves subset groups alongside All merged", {
  resolved <- RiboCrypt:::observatory_resolve_library_selections(
    library_selections = list("1" = NULL, "2" = c("SRR9", "SRR10")),
    all_library_runs = c("SRR1", "SRR2", "SRR9", "SRR10"),
    active_selection_id = "2",
    library_selection_labels = list("1" = "All merged", "2" = "subset")
  )

  expect_equal(resolved[["1"]], c("SRR1", "SRR2", "SRR9", "SRR10"))
  expect_equal(resolved[["2"]], c("SRR9", "SRR10"))
})

test_that("observatory browser restores missing All merged group from selection index", {
  resolved <- RiboCrypt:::observatory_resolve_library_selections(
    library_selections = list("2" = c("SRR9", "SRR10")),
    all_library_runs = c("SRR1", "SRR2", "SRR9", "SRR10"),
    active_selection_id = "2",
    library_selection_labels = list("1" = "All merged", "2" = "subset"),
    selection_index = c("1", "2")
  )

  expect_equal(resolved[["1"]], c("SRR1", "SRR2", "SRR9", "SRR10"))
  expect_equal(resolved[["2"]], c("SRR9", "SRR10"))
})

test_that("browser collection controller keeps All merged alongside subset group", {
  tx <- GenomicRanges::GRangesList(
    tx1 = GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 10), "+")
  )
  display_region <- tx
  experiment_df <- structure(data.frame(dummy = 1), experiment = "exp-a")
  reads_dt <- data.table::data.table(
    SRR1 = c(1, 2, 3),
    SRR2 = c(4, 5, 6),
    SRR9 = c(7, 8, 9),
    SRR10 = c(10, 11, 12)
  )
  tmp_dir <- withr::local_tempdir()
  forward_fst <- file.path(tmp_dir, "dummy-forward.fst")
  index_fst <- file.path(tmp_dir, "dummy-index.fst")
  fst::write_fst(as.data.frame(reads_dt), forward_fst)
  fst::write_fst(data.frame(file_forward = forward_fst, dummy = 1), index_fst)

  local_mocked_bindings(
    collection_path_from_exp = function(experiment_df, selected_tx, grl_all = tx()) index_fst,
    load_collection = function(path, grl = NULL, columns = NULL) {
      as.data.frame(reads_dt[, ..columns])
    },
    smoothenMultiSampCoverage = function(x, kmer, kmers_type = "mean", split_by_frame = TRUE) x,
    .package = "RiboCrypt"
  )

  resolved <- RiboCrypt:::observatory_resolve_library_selections(
    library_selections = list("1" = NULL, "2" = c("SRR9", "SRR10")),
    all_library_runs = c("SRR1", "SRR2", "SRR9", "SRR10"),
    active_selection_id = "2",
    library_selection_labels = list("1" = "All merged", "2" = "subset")
  )

  out <- RiboCrypt:::browser_collection_controller_data(
    input = list(
      extendLeaders = 0,
      extendTrailers = 0,
      kmer = 1
    ),
    selected_tx = "tx1",
    display_region = display_region,
    experiment_df = experiment_df,
    tx = function() tx,
    library_selections = resolved,
    library_selection_labels = list("1" = "All merged", "2" = "subset")
  )

  expect_equal(length(out$library_selections), 2)
  expect_true("All merged" %in% names(out$library_selections))
  expect_true(any(grepl("subset", names(out$library_selections), fixed = TRUE)))
  expect_equal(sort(out$library_selections[["All merged"]]), c("SRR1", "SRR10", "SRR2", "SRR9"))
  expect_equal(sort(out$library_selections[[grep("subset", names(out$library_selections), value = TRUE)[1]]]), c("SRR10", "SRR9"))
})
