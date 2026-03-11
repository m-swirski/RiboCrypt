make_profile_observatory_selector_fixture <- function(n = 5000L) {
  runs <- sprintf("SRR%05d", seq_len(n))
  tissues <- c("brain", "heart", "liver", "lung")
  cell_lines <- c("A1", "A2", "H1", "L1", "Lu1")
  projects <- sprintf("PRJ%03d", rep(seq_len(max(2L, ceiling(n / 50L))), length.out = n))
  authors <- sprintf("Author%02d", rep(seq_len(20L), length.out = n))
  inhibitors <- c("none", "drugA", "drugB")

  libraries_df <- data.table::data.table(
    Run = runs,
    LIBRARYTYPE = "RFP",
    tissue = rep(tissues, length.out = n),
    cell_line = rep(cell_lines, length.out = n),
    BioProject = projects,
    author = authors,
    inhibitors = rep(inhibitors, length.out = n)
  )
  libraries_df[, tissues_cell_lines := paste(tissue, cell_line, sep = " | ")]

  umap_df <- data.table::copy(libraries_df)[, `:=`(
    sample = Run,
    color_column = tissue,
    `UMAP 1` = seq_len(n) %% 97,
    `UMAP 2` = seq_len(n) %% 89
  )]
  data.table::setattr(umap_df, "color.by", "tissue")

  list(
    all_exp = data.frame(name = c("exp-a", "exp-b"), stringsAsFactors = FALSE),
    libraries_df = libraries_df,
    umap_df = umap_df,
    experiment_df = ORFik::ORFik.template.experiment()[9:10, ],
    browser_options = c(default_experiment_meta = "exp-a")
  )
}

profile_observatory_selector_setup <- function(
  n = 5000L,
  out_html = tempfile("observatory-selector-profile-", fileext = ".html")
) {
  fixture <- make_profile_observatory_selector_fixture(n)

  harness <- function(id, all_exp, experiment_df, libraries_df, browser_options = c()) {
    shiny::moduleServer(id, function(input, output, session) {
      selections <- RiboCrypt:::observatory_selector_server(
        "selector",
        all_exp,
        shiny::reactive(experiment_df),
        libraries_df,
        experiments = NULL,
        org = shiny::reactive(NULL),
        rv = shiny::reactiveValues(),
        browser_options = browser_options,
        observatory_url_state = shiny::reactiveVal(NULL)
      )

      selected_libraries <- selections$selected_libraries
      active_selection_id <- selections$active_selection_id
    })
  }

  prof <- profvis::profvis({
    testthat::local_mocked_bindings(
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
      harness,
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
        session$setInputs(
          `selector-libraries_umap_plot_selection` = fixture$libraries_df$Run[1:min(250L, n)]
        )
        session$flushReact()
        invisible(list(
          active_selection_id = active_selection_id(),
          selection_size = length(selected_libraries$data_table_selections()[["1"]])
        ))
      }
    )
  })

  htmlwidgets::saveWidget(prof, out_html, selfcontained = FALSE)
  out_html
}
