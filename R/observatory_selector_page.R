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
        3,
        library_selection_picker(ns("library_selection")),
        offset = 4
      )#,
      # shiny::column(1, library_selection_reset_button(ns("library_selection")))
    ),
    shiny::fluidRow(
      shiny::fluidRow(
        plotly::plotlyOutput(ns("libraries_umap_plot"), height = "600px") |>
          shinycssloaders::withSpinner(color = "#0dc5c1")
      ),
      shiny::fluidRow(
        DT::DTOutput(ns("libraries_data_table")) |>
          shinycssloaders::withSpinner(color = "#0dc5c1")
      )
    )
  )
}

observatory_selection_subset_label <- function(selected_runs, libraries_df) {
  if (is.null(selected_runs) || length(selected_runs) == 0 || is.null(libraries_df) || nrow(libraries_df) == 0) {
    return("")
  }

  selected_dt <- libraries_df[Run %in% as.character(selected_runs)]
  if (nrow(selected_dt) == 0) return("")

  first_unique_value <- function(column) {
    if (!(column %in% colnames(selected_dt))) return(NULL)
    values <- unique(as.character(selected_dt[[column]]))
    values <- values[!is.na(values) & nzchar(values)]
    if (length(values) == 1) values[[1]] else NULL
  }

  label_base <- first_unique_value("BioProject")
  if (is.null(label_base)) label_base <- first_unique_value("CELL_LINE")
  if (is.null(label_base)) label_base <- first_unique_value("TISSUE")
  if (is.null(label_base)) return("")

  paste(label_base, "subset")
}

observatory_selector_server <- function(
  id,
  meta_experiment_list,
  meta_experiment_df,
  all_libraries_df, experiments, org, rv, browser_options,
  observatory_url_state = shiny::reactiveVal(NULL)
) {
  shiny::moduleServer(id, function(input, output, session) {
    allsamples_observer_controller(input, output, session)

    obs_state <- shiny::reactive({
      observatory_url_state()
    }) |> shiny::bindEvent(observatory_url_state())

    url_state_pending_autogo <- shiny::reactiveVal(NULL)
    url_state_autogo_counter <- shiny::reactiveVal(0L)

    shiny::observe({
      st <- obs_state()
      if (is.null(st)) return()
      if (!is.null(st$exp) && st$exp %in% meta_experiment_list$name) {
        shiny::updateSelectizeInput(session, "dff", selected = st$exp)
      }
      if (!is.null(st$color_by) && length(st$color_by) > 0) {
        shiny::updateSelectInput(session, "color_by", selected = st$color_by)
      }
      url_state_pending_autogo(st)
    }) |> shiny::bindEvent(obs_state(), ignoreInit = FALSE, once = TRUE)

    shiny::observe({
      st <- url_state_pending_autogo()
      if (is.null(st)) return()

      exp_ready <- is.null(st$exp) || identical(as.character(input$dff), as.character(st$exp))
      color_ready <- is.null(st$color_by) || setequal(as.character(input$color_by), as.character(st$color_by))
      if (!exp_ready || !color_ready) return()

      url_state_autogo_counter(url_state_autogo_counter() + 1L)
      url_state_pending_autogo(NULL)
    }) |> shiny::bindEvent(
      url_state_pending_autogo(),
      input$dff,
      input$color_by,
      ignoreInit = TRUE
    )

    observatory_module <- shiny::reactive({
      create_observatory_module(meta_experiment_df(), all_libraries_df)
    }) |> shiny::bindCache(name(meta_experiment_df())) |> shiny::bindEvent(
      input$go,
      url_state_autogo_counter()
    )

    libraries_df <- shiny::reactive({
      observatory_module()$get_libraries_data()
    })

    output$libraries_umap_plot <- plotly::renderPlotly({
      time_before <- Sys.time()
      p <- observatory_module()$get_umap_data(input$color_by) |>
        umap_plot(color.by = input$color_by) |>
        htmlwidgets::onRender(
          fetchJS("umap_plot_extension.js"),
          session$ns("libraries_umap_plot_selection")
        )
      timer_done_nice_print("-- UMAP obs plotly done: ", time_before)
      p
    }) |> shiny::bindCache(name(meta_experiment_df()), input$color_by) |> shiny::bindEvent(observatory_module())

    plot_selection <- shiny::reactive({
      shiny::req(!is.null(input$libraries_umap_plot_selection))
      input$libraries_umap_plot_selection
    })

    output$libraries_data_table <- DT::renderDT({
      DT::datatable(
        filtered_libraries_df(),
        filter = "top",
        options = list(
          searching = TRUE,
          pageLength = 25,
          initComplete = DT::JS(
            "function(settings, json) {",
            "  var table = this.api();",
            "  var $container = $(table.table().container());",
            "  var baseId = table.table().node().id || 'libraries_data_table';",
            "  var selectedRuns = new Set();",
            "  var syncingSelection = false;",
            "  function runColumnIndex() {",
            "    var columns = table.settings()[0].aoColumns || [];",
            "    for (var i = 0; i < columns.length; i++) {",
            "      if (columns[i] && columns[i].sTitle === 'Run') return i;",
            "    }",
            "    return 0;",
            "  }",
            "  function sendSelectedRuns() {",
            "    Shiny.setInputValue(baseId + '_selected_runs', Array.from(selectedRuns), {priority: 'event'});",
            "  }",
            "  function sendFilters() {",
            "    var global = $container.find('div.dataTables_filter input').val() || '';",
            "    var columnFilters = [];",
            "    $container.find('thead tr input').each(function(){",
            "      columnFilters.push($(this).val() || '');",
            "    });",
            "    Shiny.setInputValue(baseId + '_manual_search', global, {priority: 'event'});",
            "    Shiny.setInputValue(baseId + '_manual_search_columns', columnFilters, {priority: 'event'});",
            "  }",
            "  function applySelectionToCurrentPage() {",
            "    var runCol = runColumnIndex();",
            "    syncingSelection = true;",
            "    table.rows({page: 'current'}).every(function() {",
            "      var data = this.data();",
            "      var run = Array.isArray(data) ? data[runCol] : data.Run;",
            "      if (selectedRuns.has(run)) {",
            "        this.select();",
            "      } else {",
            "        this.deselect();",
            "      }",
            "    });",
            "    syncingSelection = false;",
            "  }",
            "  function syncSetFromSelection(type, indexes) {",
            "    if (syncingSelection) return;",
            "    var runCol = runColumnIndex();",
            "    (indexes || []).forEach(function(idx) {",
            "      var data = table.row(idx).data();",
            "      var run = Array.isArray(data) ? data[runCol] : data.Run;",
            "      if (!run) return;",
            "      if (type === 'select') selectedRuns.add(run);",
            "      if (type === 'deselect') selectedRuns.delete(run);",
            "    });",
            "    sendSelectedRuns();",
            "  }",
            "  $container.on('keyup change', 'div.dataTables_filter input', sendFilters);",
            "  $container.on('keyup change', 'thead tr input', sendFilters);",
            "  table.on('select.dt', function(e, dt, type, indexes) {",
            "    if (type === 'row') syncSetFromSelection('select', indexes);",
            "  });",
            "  table.on('deselect.dt', function(e, dt, type, indexes) {",
            "    if (type === 'row') syncSetFromSelection('deselect', indexes);",
            "  });",
            "  table.on('draw.dt', function() {",
            "    applySelectionToCurrentPage();",
            "    sendFilters();",
            "  });",
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
            "  Shiny.addCustomMessageHandler('librariesApplyTableSelection', function(message) {",
            "    if (!message || message.table_id !== baseId) return;",
            "    selectedRuns = new Set(message.selected_runs || []);",
            "    applySelectionToCurrentPage();",
            "    sendSelectedRuns();",
            "  });",
            "}"
          )
        ),
        selection = "multiple"
      )
    }, server = TRUE)

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
      data_table_selection_val,
      initial_state = shiny::reactive({
        st <- obs_state()
        if (is.null(st) || is.null(st$selections)) return(NULL)
        st$selections
      })
    )

    filtered_libraries_df <- shiny::reactive({
      libs <- libraries_df()
      plot_sel <- input$libraries_umap_plot_selection
      if (is.null(plot_sel) || length(plot_sel) == 0) {
        return(libs)
      }
      libs[Run %in% plot_sel]
    }) |> shiny::bindEvent(
      libraries_df(),
      input$libraries_umap_plot_selection,
      ignoreNULL = FALSE
    )

    current_table_filters <- shiny::reactive({
      manual_global <- input$libraries_data_table_manual_search
      manual_columns <- input$libraries_data_table_manual_search_columns
      if (is.null(manual_global) && is.null(manual_columns)) {
        manual_global <- input$libraries_data_table_search
        manual_columns <- input$libraries_data_table_search_columns
      }
      list(
        global = manual_global %||% "",
        columns = manual_columns %||% rep("", ncol(libraries_df()))
      )
    })

    data_table_filtered_df <- shiny::reactive({
      filters <- current_table_filters()
      apply_dt_filters(
        filtered_libraries_df(),
        filters$global,
        filters$columns
      )
    })

    data_table_selection <- shiny::reactive({
      filtered_df <- data_table_filtered_df()
      selected_runs <- input$libraries_data_table_selected_runs

      if (is.null(selected_runs) || length(selected_runs) == 0) {
        return(filtered_df$Run)
      }

      intersect(filtered_df$Run, as.character(selected_runs))
    })

    shiny::observeEvent(
      current_table_filters(),
      {
        selection_value <- tryCatch(data_table_selection(), error = function(e) NULL)
        if (!is.null(selection_value)) {
          data_table_selection_val(selection_value)
        }
      },
      ignoreInit = FALSE
    )

    shiny::observeEvent(
      filtered_libraries_df(),
      {
        selection_value <- tryCatch(data_table_selection(), error = function(e) NULL)
        if (!is.null(selection_value)) {
          data_table_selection_val(selection_value)
        }
      },
      ignoreInit = FALSE
    )

    shiny::observeEvent(
      input$libraries_data_table_selected_runs,
      {
        selection_value <- tryCatch(data_table_selection(), error = function(e) NULL)
        if (!is.null(selection_value)) {
          data_table_selection_val(selection_value)
        }
      },
      ignoreInit = FALSE
    )

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
          table_id = session$ns("libraries_data_table"),
          global = filters$global %||% "",
          columns = filters$columns
        )
      )
      invisible(NULL)
    }

    apply_data_table_selection <- function(selection_id) {
      selected_runs <- selected_libraries$all_selections$data_table_selections()[[selection_id]]
      if (is.null(selected_runs)) selected_runs <- character()
      session$sendCustomMessage(
        "librariesApplyTableSelection",
        list(
          table_id = session$ns("libraries_data_table"),
          selected_runs = as.character(selected_runs)
        )
      )
      invisible(NULL)
    }

    shiny::observe({
      shiny::req(selected_libraries$active_selection_id())
      plot_sel <- input$libraries_umap_plot_selection
      if (is.null(plot_sel) || length(plot_sel) != 0) {
        return()
      }

      selection_id <- selected_libraries$active_selection_id()
      filters <- data_table_filters()
      filters[[selection_id]] <- list(
        global = "",
        columns = rep("", ncol(libraries_df()))
      )
      data_table_filters(filters)
      DT::updateSearch(
        libraries_data_table_proxy(),
        keywords = list(
          global = "",
          columns = rep("", ncol(libraries_df()))
        )
      )
      apply_data_table_filters(selection_id)
    }) |> shiny::bindEvent(
      input$libraries_umap_plot_selection,
      ignoreInit = TRUE,
      ignoreNULL = FALSE
    )

    shiny::observe({
      shiny::req(selected_libraries$active_selection_id())
      plot_sel <- input$libraries_umap_plot_selection
      if (is.null(plot_sel) || length(plot_sel) == 0) return()

      selected_libraries$set_active_label(
        observatory_selection_subset_label(plot_sel, libraries_df())
      )
    }) |> shiny::bindEvent(
      input$libraries_umap_plot_selection,
      ignoreInit = TRUE,
      ignoreNULL = FALSE
    )

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
      apply_data_table_selection(selection_id)
    }) |> shiny::bindEvent(selected_libraries$active_selection_id())

    shiny::observe({
      columns <- input$libraries_data_table_manual_search_columns
      global <- input$libraries_data_table_manual_search
      if (is.null(columns) && is.null(global)) {
        columns <- input$libraries_data_table_search_columns
        global <- input$libraries_data_table_search
      }
      if (is.null(columns) && is.null(global)) {
        selected_libraries$set_active_label("")
        return()
      }

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
      if (length(terms) == 0) {
        selected_libraries$set_active_label("")
        return()
      }

      label <- paste(tolower(unique(terms)), collapse = "|")
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
      shiny::req(selected_libraries$active_selection_id())
      filters <- current_table_filters()
      DT::replaceData(
        libraries_data_table_proxy(),
        filtered_libraries_df(),
        resetPaging = FALSE,
        clearSelection = "none"
      )
      DT::updateSearch(
        libraries_data_table_proxy(),
        keywords = list(
          global = filters$global,
          columns = filters$columns
        )
      )
      apply_data_table_selection(selected_libraries$active_selection_id())
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
      selected_libraries = selected_libraries$all_selections,
      selected_experiment = shiny::reactive(input$dff),
      color_by = shiny::reactive(input$color_by),
      active_selection_id = selected_libraries$active_selection_id
    )
  })
}
