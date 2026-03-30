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
        DT::DTOutput(ns("libraries_data_table"))
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

observatory_priority_library_columns <- function(df) {
  if (is.null(df) || ncol(df) == 0) {
    return(character())
  }

  priority_columns <- c(
    "Run",
    "BioProject",
    "TISSUE",
    "CELL_LINE",
    "tissue",
    "cell_line",
    "author",
    "inhibitors",
    "LIBRARYTYPE"
  )

  priority_columns <- intersect(priority_columns, colnames(df))
  c(priority_columns, setdiff(colnames(df), priority_columns))
}

observatory_format_libraries_df <- function(df, digits = 3L) {
  if (is.null(df) || nrow(df) == 0) {
    return(df)
  }

  formatted_df <- data.table::copy(df)
  drop_columns <- intersect(c("LIBRARYTYPE", "ScientificName"), colnames(formatted_df))
  if (length(drop_columns) > 0) {
    formatted_df[, (drop_columns) := NULL]
  }
  data.table::setcolorder(formatted_df, observatory_priority_library_columns(formatted_df))
  formatted_df
}

observatory_selection_is_all_merged <- function(selected_runs, all_runs) {
  if (is.null(all_runs)) {
    return(is.null(selected_runs) || length(selected_runs) == 0)
  }

  all_runs <- unique(as.character(all_runs))
  if (length(all_runs) == 0) {
    return(is.null(selected_runs) || length(selected_runs) == 0)
  }

  if (is.null(selected_runs)) {
    return(TRUE)
  }

  selected_runs <- unique(as.character(selected_runs))
  setequal(selected_runs, all_runs)
}

observatory_selector_umap_plot_shiny <- function(observatory_module, color_by, session) {
  time_before <- Sys.time()
  p <- observatory_module$get_umap_data(color_by) |>
    umap_plot(color.by = color_by) |>
    htmlwidgets::onRender(
      fetchJS("umap_plot_extension.js"),
      session$ns("libraries_umap_plot_selection")
    )
  timer_done_nice_print("-- UMAP obs plotly done: ", time_before)
  p
}

observatory_selector_data_table_shiny <- function(libraries_df) {
  time_before <- Sys.time()
  table_df <- observatory_format_libraries_df(libraries_df)
  numeric_targets <- unname(which(vapply(table_df, is.numeric, logical(1))) - 1L)
  left_columns <- min(sum(colnames(table_df) %in% c("Run", "BioProject", "TISSUE")), 3L)

  dt <- DT::datatable(
    table_df,
    filter = "top",
    extensions = "FixedColumns",
    options = list(
      processing = FALSE,
      searching = TRUE,
      pageLength = 15,
      lengthMenu = c(15, 25, 50, 100),
      scrollX = TRUE,
      autoWidth = FALSE,
      deferRender = TRUE,
      fixedColumns = list(leftColumns = left_columns),
      columnDefs = list(
        list(
          targets = numeric_targets,
          className = "dt-body-right dt-head-right"
        )
      ),
      initComplete = DT::JS(
        "function(settings, json) {",
        "  var table = this.api();",
        "  var $container = $(table.table().container());",
        "  var baseId = table.table().node().id || 'libraries_data_table';",
        sprintf("  var numericColumns = new Set(%s);", jsonlite::toJSON(as.integer(numeric_targets), auto_unbox = TRUE)),
        "  var selectedRuns = new Set();",
        "  var syncingSelection = false;",
        "  function normalizeFilterLayout() {",
        "    $container.css('overflow-x', 'auto');",
        "    $container.find('thead tr:eq(1) th').each(function(i) {",
        "      var $input = $(this).find('input');",
        "      if ($input.length === 0) return;",
        "      $input.attr('title', table.column(i).header().textContent || '');",
        "      $input.attr('placeholder', numericColumns.has(i) ? 'filter' : 'search');",
        "      if (numericColumns.has(i)) {",
        "        $input.addClass('observatory-dt-numeric-filter');",
        "      }",
        "    });",
        "  }",
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
        "    normalizeFilterLayout();",
        "    applySelectionToCurrentPage();",
        "    sendFilters();",
        "  });",
        "  normalizeFilterLayout();",
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
  timer_done_nice_print("-- DT obs init done: ", time_before)
  dt
}

observatory_apply_dt_filters <- function(df, global_search, column_searches) {
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

observatory_selector_additional_controller <- function(input, output, session, observatory) {
  `%||%` <- function(x, y) if (is.null(x)) y else x

  current_plot_selection <- observatory$current_plot_selection
  observatory_module <- observatory$observatory_module
  libraries_df <- observatory$libraries_df
  selected_libraries <- observatory$selected_libraries
  libraries_data_table_proxy <- observatory$libraries_data_table_proxy
  filtered_libraries_df <- observatory$filtered_libraries_df
  current_table_filters <- observatory$current_table_filters
  data_table_selection <- observatory$data_table_selection
  data_table_selection_val <- observatory$data_table_selection_val

  data_table_filters <- shiny::reactiveVal(list())
  last_active_selection_id <- shiny::reactiveVal(NULL)

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
    current_plot_selection(input$libraries_umap_plot_selection)
  }) |> shiny::bindEvent(
    input$libraries_umap_plot_selection,
    ignoreInit = FALSE,
    ignoreNULL = FALSE
  )

  shiny::observe({
    observatory_module()
    current_plot_selection(NULL)
    if (!is.null(selected_libraries$active_selection_id())) {
      selected_libraries$set_active_label("All merged")
    }
  }) |> shiny::bindEvent(observatory_module(), ignoreInit = TRUE)

  shiny::observeEvent(
    list(
      current_table_filters(),
      filtered_libraries_df(),
      input$libraries_data_table_selected_runs
    ),
    {
      selection_value <- tryCatch(data_table_selection(), error = function(e) NULL)
      data_table_selection_val(selection_value)
    },
    ignoreInit = FALSE
  )

  shiny::observe({
    shiny::req(selected_libraries$active_selection_id())
    plot_sel <- current_plot_selection()
    if (!is.null(plot_sel) && length(plot_sel) != 0) {
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
    selected_libraries$set_active_label("All merged")
    apply_data_table_filters(selection_id)
  }) |> shiny::bindEvent(
    current_plot_selection(),
    ignoreInit = TRUE,
    ignoreNULL = FALSE
  )

  shiny::observe({
    shiny::req(selected_libraries$active_selection_id())
    plot_sel <- current_plot_selection()
    if (is.null(plot_sel) || length(plot_sel) == 0) return()

    selected_libraries$set_active_label(
      observatory_selection_subset_label(plot_sel, libraries_df())
    )
  }) |> shiny::bindEvent(
    current_plot_selection(),
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
      selected_libraries$set_active_label("All merged")
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
      selected_libraries$set_active_label("All merged")
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
    updated_df <- observatory_format_libraries_df(filtered_libraries_df())
    DT::replaceData(
      libraries_data_table_proxy(),
      updated_df,
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
  }) |> shiny::bindEvent(filtered_libraries_df(), ignoreInit = TRUE)

  shiny::observe({
    selection <- selected_libraries$active_data_table_selection()
    if (observatory_selection_is_all_merged(selection, libraries_df()$Run)) {
      session$sendCustomMessage(
        "librariesActiveSelectionReset",
        ""
      )
    } else {
      session$sendCustomMessage(
        "librariesActiveSelectionChanged",
        selection
      )
    }
  }) |> shiny::bindEvent(
    selected_libraries$active_data_table_selection(),
    libraries_df(),
    ignoreNULL = FALSE
  )
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
      observatory_selector_umap_plot_shiny(
        observatory_module = observatory_module(),
        color_by = input$color_by,
        session = session
      )
    }) |> shiny::bindCache(name(meta_experiment_df()), input$color_by) |> shiny::bindEvent(observatory_module())

    current_plot_selection <- shiny::reactiveVal(NULL)

    plot_selection <- shiny::reactive({
      plot_sel <- current_plot_selection()
      if (is.null(plot_sel) || length(plot_sel) == 0) {
        return(libraries_df()$Run)
      }
      plot_sel
    })

    output$libraries_data_table <- DT::renderDT({
      observatory_selector_data_table_shiny(libraries_df())
    }, server = TRUE)

    libraries_data_table_proxy <- shiny::reactive({
      DT::dataTableProxy("libraries_data_table")
    })

    data_table_selection_val <- shiny::reactiveVal(NULL)

    selected_libraries <- library_selection_server(
      "library_selection",
      plot_selection,
      data_table_selection_val,
      default_selection = shiny::reactive(libraries_df()$Run),
      initial_state = shiny::reactive({
        st <- obs_state()
        if (is.null(st) || is.null(st$selections)) return(NULL)
        st$selections
      })
    )

    filtered_libraries_df <- shiny::reactive({
      libs <- libraries_df()
      plot_sel <- current_plot_selection()
      if (is.null(plot_sel) || length(plot_sel) == 0) {
        return(libs)
      }
      libs[Run %in% plot_sel]
    }) |> shiny::bindEvent(
      libraries_df(),
      current_plot_selection(),
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
      observatory_apply_dt_filters(
        observatory_format_libraries_df(filtered_libraries_df()),
        filters$global,
        filters$columns
      )
    })

    data_table_selection <- shiny::reactive({
      filtered_df <- data_table_filtered_df()
      selected_runs <- input$libraries_data_table_selected_runs
      filters <- current_table_filters()
      plot_sel <- current_plot_selection()
      has_filters <- nzchar(filters$global) || any(nzchar(filters$columns))
      has_plot_subset <- !is.null(plot_sel) && length(plot_sel) > 0

      if (is.null(selected_runs) || length(selected_runs) == 0) {
        if (!has_filters && !has_plot_subset) return(libraries_df()$Run)
        return(filtered_df$Run)
      }

      intersect(filtered_df$Run, as.character(selected_runs))
    })

    observatory_selector_additional_controller(
      input, output, session,
      observatory = list(
        current_plot_selection = current_plot_selection,
        observatory_module = observatory_module,
        libraries_df = libraries_df,
        selected_libraries = selected_libraries,
        libraries_data_table_proxy = libraries_data_table_proxy,
        filtered_libraries_df = filtered_libraries_df,
        current_table_filters = current_table_filters,
        data_table_selection = data_table_selection,
        data_table_selection_val = data_table_selection_val
      )
    )

    list(
      meta_experiment_df = meta_experiment_df,
      selected_libraries = selected_libraries$all_selections,
      all_library_runs = shiny::reactive(libraries_df()$Run),
      selected_experiment = shiny::reactive(input$dff),
      color_by = shiny::reactive(input$color_by),
      active_selection_id = selected_libraries$active_selection_id
    )
  })
}
