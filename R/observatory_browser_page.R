observatory_browser_ui <- function(id, browser_options) {
  ns <- shiny::NS(id)

  shiny::tabPanel(
    "Browse",
    shiny::fluidRow(
      shiny::column(
        2,
        shiny::selectizeInput(ns("gene"),
          "Gene",
          choices = list()
        )
      ),
      shiny::column(
        2,
        shiny::selectizeInput(
          ns("tx"),
          "Transcript",
          choices = list()
        )
      ),
      shiny::column(1, plot_button(ns("go"))),
      shiny::column(1, shiny::uiOutput(ns("clip"))),
      shiny::column(5,
        shiny::fluidRow(

          shiny::column(
            3,
            frame_type_select(ns, selected = browser_options["default_frame_type"])
          ),
          shiny::column(
            3,
            sliderInput(ns("kmer"), "K-mer length", min = 1, max = 20,
                        value = as.numeric(browser_options["default_kmer"]))
          ),
          shiny::column(
            3,
            shiny::numericInput(ns("extendLeaders"), "5' extension", 0)
          ),
          shiny::column(
            3,
            shiny::numericInput(ns("extendTrailers"), "3' extension", 0)
          )
        ),
        offset = 2
      )
    ),
    shiny::fluidRow(
      shiny::column(
        12,
        jqui_resizable(plotly::plotlyOutput(ns("browser_plot"), height = "700px")) |>
          shinycssloaders::withSpinner(color = "#0dc5c1")
      )
    )
  )
}

# TODO
# Add settings for kmer, leader and trailer extensions and aggregation method
observatory_browser_server <- function(
  id,
  df,
  library_selections,
  library_selection_labels,
  gene_name_list, experiments, org, gg_theme,
  rv, browser_options,
  selection_index = shiny::reactive(NULL),
  active_selection_id = shiny::reactive(NULL),
  selected_experiment = shiny::reactive(NULL),
  color_by = shiny::reactive(NULL),
  observatory_url_state = shiny::reactiveVal(NULL)
) {
  shiny::moduleServer(id, function(input, output, session) {
    # -- Gene / transcript input wiring ------------------------------------
    allsamples_observer_controller(input, output, session)
    kickoff <- shiny::reactiveVal(FALSE)

    output$clip <- shiny::renderUI({
      state <- observatory_state_from_inputs(
        selected_experiment = selected_experiment(),
        color_by = color_by(),
        view = "browser",
        browser = list(
          gene = input$gene,
          tx = input$tx,
          frames_type = input$frames_type,
          kmer = input$kmer,
          extendLeaders = input$extendLeaders,
          extendTrailers = input$extendTrailers,
          go = TRUE
        ),
        selections = list(
          index = selection_index(),
          plot_selections = library_selections(),
          data_table_selections = library_selections(),
          labels = library_selection_labels(),
          active_selection_id = active_selection_id()
        )
      )
      observatory_clipboard_url_button(state, session)
    })

    shiny::observe({
      st <- observatory_url_state()
      if (is.null(st)) return()
      br <- st$browser
      if (is.null(br)) return()

      if (!is.null(br$gene)) shiny::updateSelectizeInput(session, "gene", selected = br$gene)
      if (!is.null(br$tx)) shiny::updateSelectizeInput(session, "tx", selected = br$tx)
      if (!is.null(br$frames_type)) shiny::updateSelectizeInput(session, "frames_type", selected = br$frames_type)
      if (!is.null(br$kmer) && !is.na(as.numeric(br$kmer))) shiny::updateSliderInput(session, "kmer", value = as.numeric(br$kmer))
      if (!is.null(br$extendLeaders) && !is.na(as.numeric(br$extendLeaders))) {
        shiny::updateNumericInput(session, "extendLeaders", value = as.numeric(br$extendLeaders))
      }
      if (!is.null(br$extendTrailers) && !is.na(as.numeric(br$extendTrailers))) {
        shiny::updateNumericInput(session, "extendTrailers", value = as.numeric(br$extendTrailers))
      }
    }) |> shiny::bindEvent(observatory_url_state(), ignoreInit = FALSE, once = TRUE)

    shiny::observe({
      st <- observatory_url_state()
      if (is.null(st) || !identical(st$view, "browser")) return()
      if (!isTRUE(st$browser$go)) return()
      req(nzchar(input$gene), nzchar(input$tx))
      sel <- library_selections()
      if (is.null(sel) || !any(lengths(sel) > 0)) return()
      kickoff(TRUE)
    }) |> shiny::bindEvent(
      input$gene, input$tx,
      library_selections(),
      observatory_url_state(),
      ignoreInit = TRUE
    )
    # -- Annotation loading ------------------------------------------------
    # Load transcript and CDS annotation once per experiment.
    tx <- shiny::reactive({
      loadRegion(df(), "tx")
    }) |>
      shiny::bindCache(name(df())) |>
      shiny::bindEvent(df(), ignoreNULL = TRUE)

    cds <- shiny::reactive({
      loadRegion(df(), "cds")
    }) |>
      shiny::bindCache(name(df())) |>
      shiny::bindEvent(df(), ignoreNULL = TRUE)

    # -- Plot controller ---------------------------------------------------

    # Build a minimal set of plot controls when "go" is pressed.
    # Uses hardcoded defaults (no advanced settings panel in UI).
    main_plot_controls <- shiny::reactive({
      selected_tx <- input$tx
      aggregation_method <- rowMeans
      shiny::req(selected_tx, selected_tx != "")
      time_before <- Sys.time()
      message("- Obs browser controller")

      display_region <- observed_tx_annotation(selected_tx, tx)
      cds_annotation <- observed_cds_annotation(selected_tx, cds)
      tx_annotation <- observed_cds_annotation(selected_tx, cds)

      # Subset experiment to selected libraries (if any)
      experiment_df <- df()
      frames_type <- input$frames_type

      path <- collection_path_from_exp(experiment_df, selected_tx, grl_all = tx())
      display_region_grl <- ORFik::extendTrailers(
        ORFik::extendLeaders(attr(path, "range"), input$extendLeaders),
        input$extendTrailers
      )

      lib_sel <- library_selections()
      if (is.null(lib_sel)) lib_sel <- list()
      lib_sel <- lapply(lib_sel, function(selection) {
        selection <- as.character(selection)
        unique(selection[!is.na(selection) & nzchar(selection)])
      })
      lib_sel <- lib_sel[lengths(lib_sel) > 0]
      if (length(lib_sel) == 0) {
        stop("No non-empty library selection groups available for observatory browse.")
      }

      runs <- unique(unlist(lib_sel, use.names = FALSE))
      reads <- load_collection(path, grl = display_region_grl, columns = runs)
      available_runs <- colnames(reads)
      lib_sel <- lapply(lib_sel, intersect, y = available_runs)
      lib_sel <- lib_sel[lengths(lib_sel) > 0]
      if (length(lib_sel) == 0) {
        stop("No selected runs are available for the chosen transcript in observatory browse.")
      }
      with_frames <- ORFik::libraryTypes(
        experiment_df,
        uniqueTypes = FALSE
      ) %in%
        c("RFP", "RPF", "LSU", "TI")

      if (!is.null(lib_sel) && !is.null(library_selection_labels)) {
        labels <- library_selection_labels()
        selection_ids <- names(lib_sel)
        display_labels <- vapply(
          selection_ids,
          function(selection_id) {
            label <- labels[[selection_id]]
            if (is.null(label) || label == "" || label == selection_id) {
              selection_id
            } else {
              paste(selection_id, label, sep = " - ")
            }
          },
          character(1)
        )
        names(lib_sel) <- display_labels
      }
      fst_index <- path
      file_forward <- fst::read_fst(path)[1,]$file_forward
      file_forward <- file.path(dirname(fst_index), basename(file_forward))
      megafst_samples <- length(fst::metadata_fst(file_forward)$columnNames)
      group_is_all <- lengths(lib_sel) == megafst_samples
      if (any(group_is_all)) names(lib_sel)[group_is_all] <- "All merged"
      message("Number of runs used: ", length(unlist(lib_sel)),
              " (", paste(lengths(lib_sel), collapse = ", "), ")")

      profiles <- lapply(lib_sel, function(selection) {
        count <- aggregation_method(
          reads[, selection, with = FALSE],
          na.rm = TRUE
        )
        data.table::data.table(
          count = count,
          position = seq_along(count),
          library = as.factor(rep_len(1, length.out = length(count))),
          frame = as.factor(rep_len(1:3, length.out = length(count)))
        ) |>
          smoothenMultiSampCoverage(
            input$kmer,
            kmers_type = "mean",
            split_by_frame = TRUE
          )
      })
      timer_done_nice_print("- Obs browser controller done: ", time_before)
      shiny::reactiveValues(
        dff = experiment_df,
        display_region = display_region,
        customRegions = NULL,
        extendTrailers = input$extendTrailers,
        extendLeaders = input$extendLeaders,
        export_format = "svg",
        summary_track = FALSE,
        summary_track_type = "lines",
        viewMode = FALSE,
        collapsed_introns_width = 0,
        kmerLength = input$kmer,
        frames_type = frames_type,
        annotation = cds_annotation,
        tx_annotation = tx_annotation,
        reads = reads,
        custom_sequence = NULL,
        log_scale = FALSE,
        phyloP = FALSE,
        withFrames = TRUE,
        zoom_range = numeric(0),
        frames_subset = "all",
        mapability = FALSE,
        frame_colors = "R",
        colors = NULL,
        library_selections = lib_sel,
        gg_theme = gg_theme,
        is_cellphone = FALSE,
        user_browser_width = NULL,
        hash_bottom = paste(selected_tx, collapse = "|"),
        hash_browser = paste(
          selected_tx, frames_type,
          paste(experiment_df$Run, collapse = ","),
          collapse = "|"
        ),
        hash_expression = paste(selected_tx, collapse = "|"),
        profiles = profiles
      )
    }) |> shiny::bindEvent(list(input$go, kickoff()), ignoreNULL = TRUE, ignoreInit = TRUE)

    # -- Plot rendering ----------------------------------------------------

    bottom_panel <- shiny::reactive({
      bottom_panel_shiny(main_plot_controls)
    }) |>
      shiny::bindCache(main_plot_controls()$hash_bottom) |>
      shiny::bindEvent(main_plot_controls(), ignoreNULL = TRUE)

    browser_plot <- shiny::reactive({
      browser_track_panel_shiny(
        main_plot_controls, bottom_panel(), session,
        ylabels = names(main_plot_controls()$library_selections),
        profiles = main_plot_controls()$profiles,
        use_fst = TRUE,
        selected_libraries = main_plot_controls()$library_selections
      )
    }) |> shiny::bindEvent(bottom_panel(), ignoreNULL = TRUE)

    output$browser_plot <- plotly::renderPlotly({
      browser_plot()
    }) |> shiny::bindEvent(browser_plot(), ignoreNULL = TRUE)
  })
}
