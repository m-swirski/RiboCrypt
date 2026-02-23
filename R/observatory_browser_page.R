observatory_browser_ui <- function(id) {
  ns <- shiny::NS(id)

  shiny::tabPanel(
    "Browse",
    shiny::fluidRow(
      shiny::column(
        2,
        shiny::selectizeInput(ns("gene_input"),
          "Gene",
          choices = list()
        )
      ),
      shiny::column(
        2,
        shiny::selectizeInput(
          ns("tx_input"),
          "Transcript",
          choices = list()
        )
      ),
      shiny::column(1, plot_button(ns("go")))
    ),
    shiny::fluidRow(
      shiny::column(
        12,
        plotly::plotlyOutput(ns("browser_plot"), height = "550px") |>
          shinycssloaders::withSpinner(color = "#0dc5c1")
      )
    )
  )
}

# TODO
# Add settings for kmer, leader and trailer extensions and aggregation method
observatory_browser_server <- function(
  id,
  meta_experiment_df,
  library_selections
) {
  shiny::moduleServer(id, function(input, output, session) {
    # -- Gene / transcript input wiring ------------------------------------

    # Derive gene <-> transcript mapping table from the current experiment.
    # Returns a data.table with columns: value (tx id), label (gene symbol).
    gene_name_list <- shiny::reactive({
      get_gene_name_categories(meta_experiment_df())
    }) |> shiny::bindEvent(meta_experiment_df(), ignoreNULL = TRUE)

    # When the experiment changes, repopulate the gene dropdown with the new
    # gene symbols and pre-select the first one.
    shiny::observe({
      gene_names <- unique(gene_name_list()[, 2][[1]])
      shiny::updateSelectizeInput(
        session,
        "gene_input",
        choices = gene_names,
        selected = gene_names[1],
        server = TRUE
      )
    }) |> shiny::bindEvent(gene_name_list(), ignoreNULL = TRUE)

    # When the selected gene changes, repopulate the transcript dropdown with
    # only the isoforms that belong to the chosen gene.
    shiny::observe({
      shiny::req(input$gene_input != "")
      isoforms <- gene_name_list()[label == input$gene_input, value]
      shiny::updateSelectizeInput(
        session,
        "tx_input",
        choices = isoforms,
        selected = isoforms[1],
        server = TRUE
      )
    }) |> shiny::bindEvent(
      input$gene_input,
      ignoreNULL = TRUE, ignoreInit = TRUE
    )

    # -- Annotation loading ------------------------------------------------

    # Load transcript and CDS annotation once per experiment.
    tx <- shiny::reactive({
      loadRegion(meta_experiment_df(), "tx")
    }) |> shiny::bindEvent(meta_experiment_df(), ignoreNULL = TRUE)

    cds <- shiny::reactive({
      loadRegion(meta_experiment_df(), "cds")
    }) |> shiny::bindEvent(meta_experiment_df(), ignoreNULL = TRUE)

    # -- Plot controller ---------------------------------------------------

    # Build a minimal set of plot controls when "go" is pressed.
    # Uses hardcoded defaults (no advanced settings panel in UI).
    main_plot_controls <- shiny::reactive({
      selected_tx <- input$tx_input
      aggregation_method <- rowMeans
      shiny::req(selected_tx, selected_tx != "")

      display_region <- observed_tx_annotation(selected_tx, tx)
      cds_annotation <- observed_cds_annotation(selected_tx, cds)
      tx_annotation <- observed_cds_annotation(selected_tx, cds)

      # Subset experiment to selected libraries (if any)
      experiment_df <- meta_experiment_df()

      path <- collection_path_from_exp(experiment_df, selected_tx)
      runs <- unlist(library_selections())
      reads <- load_collection(path, format = "wide", columns = runs)
      with_frames <- ORFik::libraryTypes(
        experiment_df,
        uniqueTypes = FALSE
      ) %in%
        c("RFP", "RPF", "LSU", "TI")

      profiles <- lapply(library_selections(), function(selection) {
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
            9,
            kmers_type = "mean",
            split_by_frame = TRUE
          )
      })


      shiny::reactiveValues(
        dff = experiment_df,
        display_region = display_region,
        customRegions = NULL,
        extendTrailers = 0,
        extendLeaders = 0,
        export_format = "svg",
        summary_track = FALSE,
        summary_track_type = "lines",
        viewMode = FALSE,
        collapsed_introns_width = 0,
        kmerLength = 9,
        frames_type = "area",
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
        gg_theme = gg_theme_template(),
        is_cellphone = FALSE,
        user_browser_width = NULL,
        hash_bottom = paste(selected_tx, collapse = "|"),
        hash_browser = paste(
          selected_tx,
          paste(experiment_df$Run, collapse = ","),
          collapse = "|"
        ),
        hash_expression = paste(selected_tx, collapse = "|"),
        profiles = profiles
      )
    }) |> shiny::bindEvent(input$go, ignoreNULL = TRUE, ignoreInit = TRUE)

    # -- Plot rendering ----------------------------------------------------

    bottom_panel <- shiny::reactive({
      bottom_panel_shiny(main_plot_controls)
    }) |> shiny::bindEvent(main_plot_controls(), ignoreNULL = TRUE)

    browser_plot <- shiny::reactive({
      browser_track_panel_shiny(
        main_plot_controls, bottom_panel(), session,
        ylabels = names(library_selections()),
        profiles = main_plot_controls()$profiles,
        use_fst = TRUE,
        selected_libraries = library_selections()
      )
    }) |> shiny::bindEvent(bottom_panel(), ignoreNULL = TRUE)

    output$browser_plot <- plotly::renderPlotly({
      browser_plot()
    }) |> shiny::bindEvent(browser_plot(), ignoreNULL = TRUE)
  })
}
