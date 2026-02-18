observatory_ui <- function(id, meta_experiment_list) {
  umap_column_names <- c(
    "Tissue" = "tissue",
    "Cell line" = "cell_line",
    "Inhibitor" = "inhibitors",
    "BioProject" = "BioProject",
    "Author" = "author"
  )
  ns <- shiny::NS(id)

  shiny::tabPanel(
    "Observatory",
    shiny::fluidRow(
      shiny::column(
        2,
        shiny::selectInput(
          ns("meta_experiment"),
          "Organism",
          meta_experiment_list$name
        )
      ),
      shiny::column(
        2,
        shiny::selectInput(
          ns("color_by"),
          "Color by",
          umap_column_names,
          selected = umap_column_names[1:2],
          multiple = TRUE
        )
      ),
      shiny::column(1, plot_button(ns("go"))),
      shiny::column(
        2,
        sample_selection_picker(ns("sample_selection")),
        offset = 4
      ),  
      shiny::column(1, sample_selection_reset_button(ns("sample_selection")))
    ),
    shiny::fluidRow(
      shiny::fluidRow(
        plotly::plotlyOutput(ns("samples_umap_plot")) |>
          shinycssloaders::withSpinner(color = "#0dc5c1")
      ),
      shiny::fluidRow(
        DT::DTOutput(ns("samples_data_table")) |>
          shinycssloaders::withSpinner(color = "#0dc5c1")
      )
    )
  )
}

observatory_server <- function(id, meta_experiment_list, all_samples_df) {
  shiny::moduleServer(id, function(input, output, session) {
    observatory_module <- shiny::reactive({
      meta_experiment_df <- ORFik::read.experiment(
        meta_experiment_list[name == input$meta_experiment][[1]],
        validate = FALSE
      )
      create_observatory_module(meta_experiment_df, all_samples_df)
    }) |> shiny::bindEvent(input$go)

    samples_df <- shiny::reactive({
      observatory_module()$get_samples_data()
    })

    output$samples_umap_plot <- plotly::renderPlotly({
      observatory_module()$get_umap_data(input$color_by) |>
        umap_plot(color.by = input$color_by) |>
        htmlwidgets::onRender(
          fetchJS("umap_plot_extension.js"),
          session$ns("samples_umap_plot_selection")
        )
    }) |> shiny::bindEvent(observatory_module())

    plot_selection <- shiny::reactive({
      shiny::req(!is.null(input$samples_umap_plot_selection))
      input$samples_umap_plot_selection
    })

    output$samples_data_table <- DT::renderDT({
      samples_df()
    })

    samples_data_table_proxy <- shiny::reactive({
      DT::dataTableProxy("samples_data_table")
    })

    data_table_selection <- shiny::reactive({
      selection_is_defined <-
        !is.null(input$samples_data_table_rows_selected)

      selected_indexes <- if (!selection_is_defined) {
        input$samples_data_table_rows_all
      } else {
        input$samples_data_table_rows_selected
      }

      shiny::isolate(filtered_samples_df())[selected_indexes]$Run
    })

    selected_samples <- sample_selection_server(
      "sample_selection",
      plot_selection,
      data_table_selection
    )

    filtered_samples_df <- shiny::reactive({
      if (!is.null(selected_samples$active_plot_selection())) {
        samples_df()[
          Run %in% selected_samples$active_plot_selection()
        ]
      } else {
        samples_df()
      }
    }) |> shiny::bindEvent(
      selected_samples$active_plot_selection(),
      ignoreNULL = FALSE
    )

    shiny::observe({
      DT::replaceData(
        samples_data_table_proxy(),
        filtered_samples_df(),
        resetPaging = FALSE
      )
    }) |> shiny::bindEvent(filtered_samples_df())


    # TODO
    # ignoreNULL false
    # switch on NULL value and perform clean up on js side
    shiny::observe({
      if (is.null(selected_samples$active_data_table_selection())) {
        session$sendCustomMessage(
          "samplesActiveSelectionReset",
          ""
        )
      } else {
        session$sendCustomMessage(
          "samplesActiveSelectionChanged",
          selected_samples$active_data_table_selection()
        )
      }
    }) |> shiny::bindEvent(selected_samples$active_data_table_selection(), ignoreNULL = FALSE)
  })
}
