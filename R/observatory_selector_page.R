observatory_selector_ui <- function(id, meta_experiment_list) {
  umap_column_names <- c(
    "Tissue" = "tissue",
    "Cell line" = "cell_line",
    "Inhibitor" = "inhibitors",
    "BioProject" = "BioProject",
    "Author" = "author"
  )
  ns <- shiny::NS(id)

  shiny::tabPanel(
    "Select libraries",
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
        library_selection_picker(ns("library_selection")),
        offset = 4
      ),
      shiny::column(1, library_selection_reset_button(ns("library_selection")))
    ),
    shiny::fluidRow(
      shiny::fluidRow(
        plotly::plotlyOutput(ns("libraries_umap_plot")) |>
          shinycssloaders::withSpinner(color = "#0dc5c1")
      ),
      shiny::fluidRow(
        DT::DTOutput(ns("libraries_data_table")) |>
          shinycssloaders::withSpinner(color = "#0dc5c1")
      )
    )
  )
}

observatory_selector_server <- function(
  id,
  meta_experiment_list,
  all_libraries_df
) {
  shiny::moduleServer(id, function(input, output, session) {
    meta_experiment_df <- shiny::reactive({
      ORFik::read.experiment(
        meta_experiment_list[name == input$meta_experiment][[1]],
        validate = FALSE
      )
    }) |> shiny::bindEvent(input$go)

    observatory_module <- shiny::reactive({
      create_observatory_module(meta_experiment_df(), all_libraries_df)
    }) |> shiny::bindEvent(meta_experiment_df())

    libraries_df <- shiny::reactive({
      observatory_module()$get_libraries_data()
    })

    output$libraries_umap_plot <- plotly::renderPlotly({
      observatory_module()$get_umap_data(input$color_by) |>
        umap_plot(color.by = input$color_by) |>
        htmlwidgets::onRender(
          fetchJS("umap_plot_extension.js"),
          session$ns("libraries_umap_plot_selection")
        )
    }) |> shiny::bindEvent(observatory_module())

    plot_selection <- shiny::reactive({
      shiny::req(!is.null(input$libraries_umap_plot_selection))
      input$libraries_umap_plot_selection
    })

    output$libraries_data_table <- DT::renderDT({
      libraries_df()
    })

    libraries_data_table_proxy <- shiny::reactive({
      DT::dataTableProxy("libraries_data_table")
    })

    data_table_selection <- shiny::reactive({
      selection_is_defined <-
        !is.null(input$libraries_data_table_rows_selected)

      selected_indexes <- if (!selection_is_defined) {
        input$libraries_data_table_rows_all
      } else {
        input$libraries_data_table_rows_selected
      }

      shiny::isolate(filtered_libraries_df())[selected_indexes]$Run
    })

    selected_libraries <- library_selection_server(
      "library_selection",
      plot_selection,
      data_table_selection
    )

    filtered_libraries_df <- shiny::reactive({
      if (!is.null(selected_libraries$active_plot_selection())) {
        libraries_df()[
          Run %in% selected_libraries$active_plot_selection()
        ]
      } else {
        libraries_df()
      }
    }) |> shiny::bindEvent(
      selected_libraries$active_plot_selection(),
      ignoreNULL = FALSE
    )

    shiny::observe({
      DT::replaceData(
        libraries_data_table_proxy(),
        filtered_libraries_df(),
        resetPaging = FALSE
      )
    }) |> shiny::bindEvent(filtered_libraries_df())

    shiny::observe({
      if (is.null(selected_libraries$active_data_table_selection())) {
        session$sendCustomMessage(
          "librariesActiveSelectionReset",
          ""
        )
      } else {
        session$sendCustomMessage(
          "librariesActiveSelectionChanged",
          selected_libraries$active_data_table_selection()
        )
      }
    }) |> shiny::bindEvent(
      selected_libraries$active_data_table_selection(),
      ignoreNULL = FALSE
    )

    list(
      meta_experiment_df = meta_experiment_df,
      selected_libraries = selected_libraries$all_selections
    )
  })
}
