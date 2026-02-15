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
          multiple = TRUE
        )
      ),
      shiny::column(1, plot_button(ns("go"))),
      shiny::column(1, plot_button(ns("go_proxy")))
    ),
    shiny::fluidRow(
      shiny::fluidRow(
        plotly::plotlyOutput(ns("samples_umap_plot")) |>
          shinycssloaders::withSpinner(color = "#0dc5c1")
      ),
      shiny::fluidRow(
        DT::DTOutput(ns("samples_datatable")) |>
          shinycssloaders::withSpinner(color = "#0dc5c1")
      )
    )
  )
}

observatory_server <- function(id, meta_experiment_list, samples_df) {
  shiny::moduleServer(id, function(input, output, session) {
    observatory_module <- shiny::reactive({
      meta_experiment_df <- ORFik::read.experiment(
        meta_experiment_list[name == input$meta_experiment][[1]],
        validate = FALSE
      )
      create_observatory_module(meta_experiment_df, samples_df)
    }) |> shiny::bindEvent(input$go)

    unfiltered_samples_df <- shiny::reactive({
      observatory_module()$get_samples_data()
    })

    filtered_samples_df <- shiny::reactiveVal()

    output$samples_umap_plot <- plotly::renderPlotly({
      observatory_module()$get_umap_data(input$color_by) |>
        umap_plot(color.by = input$color_by) |>
        htmlwidgets::onRender(fetchJS("umap_plot_extension.js"))
    }) |> shiny::bindEvent(observatory_module())

    output$samples_datatable <- DT::renderDT({
      unfiltered_samples_df()
    })

    samples_data_table_proxy <- shiny::reactive({
      DT::dataTableProxy("samples_datatable")
    })

    shiny::observe({
      DT::replaceData(
        samples_data_table_proxy(),
        filtered_samples_df(),
        resetPaging = FALSE
      )
    }) |> shiny::bindEvent(filtered_samples_df())

    datatable_selection <- shiny::reactive({
      selected_indexes <- if (is.null(input$samples_datatable_rows_selected)) {
        input$samples_datatable_rows_all
      } else {
        input$samples_datatable_rows_selected
      }

      unfiltered_samples_df()[selected_indexes]
    })
  })
}
