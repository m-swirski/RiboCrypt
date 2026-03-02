observatory_selector_ui <- function(id, meta_experiment_list, browser_options) {
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
        experiment_input_select(meta_experiment_list$name, ns, browser_options,
                                "default_experiment_meta", label = "Organism")
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
  meta_experiment_df,
  all_libraries_df, experiments, org, rv, browser_options
) {
  shiny::moduleServer(id, function(input, output, session) {
    allsamples_observer_controller(input, output, session)

    observatory_module <- shiny::reactive({
      create_observatory_module(meta_experiment_df(), all_libraries_df)
    }) |> shiny::bindCache(name(meta_experiment_df()))|> shiny::bindEvent(input$go)

    libraries_df <- shiny::reactive({
      observatory_module()$get_libraries_data()
    })

    output$libraries_umap_plot <- plotly::renderPlotly({
      observatory_module()$get_umap_data(input$color_by) |>
        umap_plot(color.by = input$color_by) |>
        plotly::layout(dragmode = "select") |>
        plotly::config(
          displaylogo = FALSE,
          modeBarButtons = list(list("select2d", "lasso2d"))
        ) |>
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
      DT::datatable(
        libraries_df(),
        filter = "top",
        options = list(
          searching = TRUE,
          initComplete = DT::JS(
            "function(settings, json) {",
            "  var table = this.api();",
            "  var $container = $(table.table().container());",
            "  var baseId = table.table().node().id || 'libraries_data_table';",
            "  function sendFilters() {",
            "    var global = $container.find('div.dataTables_filter input').val() || '';",
            "    var columnFilters = [];",
            "    $container.find('thead tr input').each(function(){",
            "      columnFilters.push($(this).val() || '');",
            "    });",
            "    Shiny.setInputValue(baseId + '_manual_search', global, {priority: 'event'});",
            "    Shiny.setInputValue(baseId + '_manual_search_columns', columnFilters, {priority: 'event'});",
            "  }",
            "  $container.on('keyup change', 'div.dataTables_filter input', sendFilters);",
            "  $container.on('keyup change', 'thead tr input', sendFilters);",
            "  sendFilters();",
            "  Shiny.addCustomMessageHandler('librariesApplyTableFilters', function(message) {",
            "    if (!message || message.table_id !== baseId) return;",
            "    var global = message.global || '';",
            "    var columns = message.columns || [];",
            "    $container.find('div.dataTables_filter input').val(global);",
            "    var $inputs = $container.find('thead tr input');",
            "    $inputs.each(function(i){ $(this).val(columns[i] || ''); });",
            "    table.search(global);",
            "    columns.forEach(function(value, idx){ table.column(idx).search(value || ''); });",
            "    table.draw();",
            "    sendFilters();",
            "  });",
            "}"
          )
        ),
        selection = "multiple"
      )
    })

    libraries_data_table_proxy <- shiny::reactive({
      DT::dataTableProxy("libraries_data_table")
    })

    apply_dt_filters <- function(df, global_search, column_searches) {
      if (is.null(df) || nrow(df) == 0) {
        return(df)
      }

      keep <- rep(TRUE, nrow(df))

      if (!is.null(column_searches)) {
        max_cols <- min(ncol(df), length(column_searches))
        for (col_index in seq_len(max_cols)) {
          pattern <- column_searches[[col_index]]
          if (!is.null(pattern) && nzchar(pattern)) {
            keep <- keep & grepl(
              pattern,
              as.character(df[[col_index]]),
              ignore.case = TRUE
            )
          }
        }
      }

      if (!is.null(global_search) && nzchar(global_search)) {
        row_match <- Reduce(
          `|`,
          lapply(df, function(col) {
            grepl(
              global_search,
              as.character(col),
              ignore.case = TRUE
            )
          })
        )
        keep <- keep & row_match
      }

      df[keep]
    }

    data_table_selection_val <- shiny::reactiveVal(NULL)

    selected_libraries <- library_selection_server(
      "library_selection",
      plot_selection,
      data_table_selection_val
    )

    filtered_libraries_df <- shiny::reactive({
      libs <- libraries_df()
      plot_sel <- input$libraries_umap_plot_selection
      message("Plot selection size: ", length(plot_sel))
      if (is.null(plot_sel) || length(plot_sel) == 0) {
        return(libs)
      }
      libs[Run %in% plot_sel]
    }) |> shiny::bindEvent(
      libraries_df(),
      input$libraries_umap_plot_selection,
      ignoreNULL = FALSE
    )

    data_table_selection <- shiny::reactive({
      message("DT selection debug: entered data_table_selection()")
      base_df <- filtered_libraries_df()
      manual_global <- input$libraries_data_table_manual_search
      manual_columns <- input$libraries_data_table_manual_search_columns
      if (is.null(manual_global) && is.null(manual_columns)) {
        manual_global <- input$libraries_data_table_search
        manual_columns <- input$libraries_data_table_search_columns
      }
      filtered_df <- apply_dt_filters(
        base_df,
        manual_global,
        manual_columns
      )
      message(
        sprintf(
          "DT filter debug: base rows=%s, filtered rows=%s, global_search='%s'",
          nrow(base_df),
          nrow(filtered_df),
          if (is.null(manual_global)) "" else manual_global
        )
      )
      selection_is_defined <-
        !is.null(input$libraries_data_table_rows_selected)

      selected_indexes <- if (!selection_is_defined) {
        seq_len(nrow(filtered_df))
      } else {
        input$libraries_data_table_rows_selected
      }

      runs <- filtered_df[selected_indexes]$Run
      message(
        sprintf(
          "DT selection debug: selected rows=%s, runs_count=%s",
          length(selected_indexes),
          length(runs)
        )
      )
      if (length(runs) > 0) {
        message(
          sprintf(
            "DT selection debug: runs (first 20)=%s",
            paste(head(runs, 20), collapse = ", ")
          )
        )
      }
      runs
    })

    shiny::observe({
      message("DT debug: observer tick (no bindEvent)")
      message(
        sprintf(
          "DT debug: search='%s', search_columns=%s",
          if (is.null(input$libraries_data_table_search)) "" else input$libraries_data_table_search,
          if (is.null(input$libraries_data_table_search_columns)) {
            "NULL"
          } else {
            paste(input$libraries_data_table_search_columns, collapse = " | ")
          }
        )
      )
      message(
        sprintf(
          "DT debug: manual_search='%s', manual_search_columns=%s",
          if (is.null(input$libraries_data_table_manual_search)) "" else input$libraries_data_table_manual_search,
          if (is.null(input$libraries_data_table_manual_search_columns)) {
            "NULL"
          } else {
            paste(input$libraries_data_table_manual_search_columns, collapse = " | ")
          }
        )
      )
      message(
        sprintf(
          "DT debug: rows_all=%s, rows_filtered=%s, rows_selected=%s",
          if (is.null(input$libraries_data_table_rows_all)) "NULL" else length(input$libraries_data_table_rows_all),
          if (is.null(input$libraries_data_table_rows_filtered)) "NULL" else length(input$libraries_data_table_rows_filtered),
          if (is.null(input$libraries_data_table_rows_selected)) "NULL" else length(input$libraries_data_table_rows_selected)
        )
      )
      selection_value <- tryCatch(data_table_selection(), error = function(e) e)
      if (inherits(selection_value, "error")) {
        message(sprintf("DT debug: data_table_selection error: %s", conditionMessage(selection_value)))
      } else {
        message(sprintf("DT debug: data_table_selection runs=%s", length(selection_value)))
        data_table_selection_val(selection_value)
      }
    })

    data_table_filters <- shiny::reactiveVal(list())
    last_active_selection_id <- shiny::reactiveVal(NULL)

    `%||%` <- function(x, y) if (is.null(x)) y else x

    apply_data_table_filters <- function(selection_id) {
      filters <- data_table_filters()[[selection_id]]
      empty_cols <- rep("", ncol(libraries_df()))
      if (is.null(filters)) {
        filters <- list(global = "", columns = empty_cols)
      } else if (is.null(filters$columns)) {
        filters$columns <- empty_cols
      }
      session$sendCustomMessage(
        "librariesApplyTableFilters",
        list(
          table_id = "libraries_data_table",
          global = filters$global %||% "",
          columns = filters$columns
        )
      )
      invisible(NULL)
    }

    shiny::observe({
      shiny::req(selected_libraries$active_selection_id())
      selection_id <- selected_libraries$active_selection_id()
      previous_id <- last_active_selection_id()
      if (!is.null(previous_id) && previous_id != selection_id) {
        filters <- data_table_filters()
        filters[[previous_id]] <- list(
          columns = input$libraries_data_table_manual_search_columns,
          global = input$libraries_data_table_manual_search
        )
        data_table_filters(filters)
      }
      last_active_selection_id(selection_id)
      apply_data_table_filters(selection_id)
    }) |> shiny::bindEvent(selected_libraries$active_selection_id())

    shiny::observe({
      columns <- input$libraries_data_table_manual_search_columns
      global <- input$libraries_data_table_manual_search
      if (is.null(columns) && is.null(global)) {
        columns <- input$libraries_data_table_search_columns
        global <- input$libraries_data_table_search
      }
      if (is.null(columns) && is.null(global)) return()

      terms <- character(0)
      if (!is.null(columns)) {
        column_values <- as.character(columns)
        terms <- unlist(strsplit(column_values, "\\|", fixed = FALSE), use.names = FALSE)
      }
      if (!is.null(global) && nzchar(global)) {
        terms <- c(terms, unlist(strsplit(global, "\\|", fixed = FALSE), use.names = FALSE))
      }

      terms <- trimws(terms)
      terms <- terms[nzchar(terms)]
      if (length(terms) == 0) return()

      label <- paste(tolower(unique(terms)), collapse = "|")
      message(sprintf("DT label debug: updating selection label to '%s'", label))
      selected_libraries$set_active_label(label)
    }) |> shiny::bindEvent(
      input$libraries_data_table_manual_search_columns,
      input$libraries_data_table_manual_search,
      input$libraries_data_table_search_columns,
      input$libraries_data_table_search
    )

    shiny::observe({
      shiny::req(selected_libraries$active_selection_id())
      selection_id <- selected_libraries$active_selection_id()
      filters <- data_table_filters()
      filters[[selection_id]] <- list(
        columns = input$libraries_data_table_manual_search_columns,
        global = input$libraries_data_table_manual_search
      )
      data_table_filters(filters)
    }) |> shiny::bindEvent(
      input$libraries_data_table_manual_search,
      input$libraries_data_table_manual_search_columns
    )

    shiny::observe({
      DT::replaceData(
        libraries_data_table_proxy(),
        filtered_libraries_df(),
        resetPaging = FALSE
      )
      apply_data_table_filters(selected_libraries$active_selection_id())
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
