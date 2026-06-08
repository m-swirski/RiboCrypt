observatory_selector_dt_style <- function() {
  shiny::tags$style(shiny::HTML("
    .observatory-dt-top {
      align-items: center;
      display: flex;
      flex-wrap: nowrap;
      font-size: 0.9rem;
      gap: 0.75rem;
      margin-bottom: 0.5rem;
      overflow-x: auto;
      white-space: nowrap;
    }

    .observatory-dt-top .dataTables_length,
    .observatory-dt-top .dataTables_filter,
    .observatory-dt-top .dataTables_info,
    .observatory-dt-top .dataTables_paginate,
    .observatory-dt-top .observatory-dt-actions,
    .observatory-dt-top .dt-length,
    .observatory-dt-top .dt-search,
    .observatory-dt-top .dt-info,
    .observatory-dt-top .dt-paging {
      clear: none !important;
      float: none !important;
      margin: 0;
      padding-top: 0;
      white-space: nowrap;
    }

    .observatory-dt-top .observatory-dt-actions {
      align-items: center;
      display: flex;
      flex: 0 0 auto;
      gap: 0.35rem;
    }

    .observatory-dt-top .observatory-dt-action {
      font-size: 0.85rem;
      line-height: 1.25;
      padding: 0.2rem 0.5rem;
    }

    .observatory-dt-top .dataTables_filter,
    .observatory-dt-top .dt-search {
      margin-left: auto;
    }

    .observatory-dt-top label {
      margin-bottom: 0;
    }

    .observatory-dt-top .pagination {
      margin: 0;
    }
  "))
}

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
    observatory_selector_dt_style(),
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

observatory_normalize_plot_selection <- function(selected_runs, all_runs) {
  if (is.null(selected_runs) || length(selected_runs) == 0) {
    return(NULL)
  }

  selected_runs <- unique(as.character(selected_runs))
  if (observatory_selection_is_all_merged(selected_runs, all_runs)) {
    return(NULL)
  }

  selected_runs
}

observatory_selector_umap_plot_shiny <- function(observatory_module, color_by, session,
                                                 templates = NULL) {
  time_before <- Sys.time()
  p <- observatory_module$get_umap_data(color_by) |>
    umap_plot(
      color.by = color_by,
      template = templates$observatory_umap_plotly
    ) |>
    htmlwidgets::onRender(
      fetchJS("umap_plot_extension.js"),
      session$ns("libraries_umap_plot_selection")
    )
  timer_done_nice_print("-- UMAP obs plotly done: ", time_before)
  p
}

observatory_selector_data_table_shiny <- function(libraries_df, table_id = "libraries_data_table") {
  time_before <- Sys.time()
  table_df <- observatory_format_libraries_df(libraries_df)
  numeric_targets <- unname(which(vapply(table_df, is.numeric, logical(1))) - 1L)
  left_columns <- min(sum(colnames(table_df) %in% c("Run", "BioProject", "TISSUE")), 3L)

  dt <- DT::datatable(
    table_df,
    rownames = FALSE,
    filter = "top",
    extensions = "FixedColumns",
    options = list(
      processing = FALSE,
      searching = TRUE,
      dom = '<"observatory-dt-top"l<"observatory-dt-actions">fip>rt',
      pageLength = 15,
      lengthMenu = c(15, 25, 50, 100),
      language = list(
        lengthMenu = "Show _MENU_ libraries",
        info = "Showing _START_ to _END_ out of _TOTAL_ libraries",
        infoEmpty = "Showing 0 to 0 out of 0 libraries",
        infoFiltered = "(filtered from _MAX_ total libraries)"
      ),
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
        sprintf("  var baseId = %s;", jsonlite::toJSON(table_id, auto_unbox = TRUE)),
        sprintf("  var numericColumns = new Set(%s);", jsonlite::toJSON(as.integer(numeric_targets), auto_unbox = TRUE)),
        "  var selectedRuns = new Set();",
        "  var syncingSelection = false;",
        "  var suppressNextFilterSend = false;",
        "  function tableHeaderInputs() {",
        "    return $(table.table().header()).find('tr:eq(1) input');",
        "  }",
        "  function tableFilterHeaderCells() {",
        "    return $(table.table().header()).find('tr:eq(1) th');",
        "  }",
        "  function globalSearchInput() {",
        "    return $container.find('div.dataTables_filter input, div.dt-search input').first();",
        "  }",
        "  function normalizeFilterLayout() {",
        "    $container.css('overflow-x', 'auto');",
        "    tableFilterHeaderCells().each(function(i) {",
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
        "    var global = globalSearchInput().val() || '';",
        "    var columnFilters = [];",
        "    tableHeaderInputs().each(function(){",
        "      columnFilters.push($(this).val() || '');",
        "    });",
        "    Shiny.setInputValue(baseId + '_manual_search', global, {priority: 'event'});",
        "    Shiny.setInputValue(baseId + '_manual_search_columns', columnFilters, {priority: 'event'});",
        "  }",
        "  function normalizeRuns(runs) {",
        "    var seen = new Set();",
        "    var out = [];",
        "    (runs || []).forEach(function(run) {",
        "      if (!run) return;",
        "      run = String(run);",
        "      if (seen.has(run)) return;",
        "      seen.add(run);",
        "      out.push(run);",
        "    });",
        "    return out;",
        "  }",
        "  function runFromData(data) {",
        "    if (!data) return null;",
        "    var runCol = runColumnIndex();",
        "    var run = Array.isArray(data) ? data[runCol] : (data.Run || data[runCol] || data[String(runCol)]);",
        "    return run == null ? null : String(run);",
        "  }",
        "  function runFromRowApi(rowApi) {",
        "    var run = runFromData(rowApi.data());",
        "    if (run) return run;",
        "    var node = rowApi.node ? rowApi.node() : null;",
        "    if (!node) return null;",
        "    run = $(node).children('td').eq(runColumnIndex()).text();",
        "    return run ? String(run) : null;",
        "  }",
        "  function runFromRowNode(row) {",
        "    var run = runFromRowApi(table.row(row));",
        "    if (run) return run;",
        "    run = $(row).children('td').eq(runColumnIndex()).text();",
        "    return run ? String(run) : null;",
        "  }",
        "  function setRowSelected(rowApi, isSelected) {",
        "    var node = rowApi && rowApi.node ? rowApi.node() : null;",
        "    if (!node) return;",
        "    $(node).toggleClass('selected', !!isSelected);",
        "  }",
        "  function collectSelectedRunsFromTable() {",
        "    var runs = Array.from(selectedRuns);",
        "    $container.find('tbody tr.selected').each(function() {",
        "      runs.push(runFromRowNode(this));",
        "    });",
        "    return normalizeRuns(runs);",
        "  }",
        "  function sendSubsetRuns(action, runs) {",
        "    Shiny.setInputValue(baseId + '_subset_runs', {",
        "      action: action,",
        "      runs: normalizeRuns(runs),",
        "      nonce: Date.now()",
        "    }, {priority: 'event'});",
        "  }",
        "  function showNoSamplesSelected() {",
        "    window.alert('No samples are selected.');",
        "  }",
        "  function setSelectedRuns(runs, action) {",
        "    runs = normalizeRuns(runs);",
        "    selectedRuns = new Set(runs);",
        "    applySelectionToCurrentPage();",
        "    if (action) sendSubsetRuns(action, runs);",
        "  }",
        "  function currentPageRuns() {",
        "    var runs = [];",
        "    table.rows({page: 'current'}).every(function() {",
        "      runs.push(runFromRowApi(this));",
        "    });",
        "    return normalizeRuns(runs);",
        "  }",
        "  function rememberClickedRowSelection(row) {",
        "    window.setTimeout(function() {",
        "      var run = runFromRowNode(row);",
        "      if (!run) return;",
        "      run = String(run);",
        "      if (selectedRuns.has(run)) selectedRuns.delete(run);",
        "      else selectedRuns.add(run);",
        "      applySelectionToCurrentPage();",
        "      sendSelectedRuns();",
        "    }, 0);",
        "  }",
        "  function installActionButtons() {",
        "    var $actions = $container.find('.observatory-dt-actions');",
        "    if ($actions.length === 0 || $actions.data('observatoryButtonsInstalled')) return;",
        "    $actions.data('observatoryButtonsInstalled', true);",
        "    var buttonClass = 'btn btn-outline-secondary btn-sm observatory-dt-action';",
        "    var $selectedOnly = $('<button/>', {type: 'button', text: 'Subset to selected', 'class': buttonClass, title: 'Use clicked table rows as this selection'});",
        "    var $pageOnly = $('<button/>', {type: 'button', text: 'Subset to page', 'class': buttonClass, title: 'Use rows on the current table page as this selection'});",
        "    var $removeSubset = $('<button/>', {type: 'button', text: 'Remove subset', 'class': buttonClass, title: 'Show all libraries again'});",
        "    $selectedOnly.on('click', function() {",
        "      var runs = collectSelectedRunsFromTable();",
        "      if (runs.length === 0) {",
        "        showNoSamplesSelected();",
        "        return;",
        "      }",
        "      setSelectedRuns(runs, 'selected');",
        "    });",
        "    $pageOnly.on('click', function() {",
        "      var runs = currentPageRuns();",
        "      if (runs.length === 0) {",
        "        showNoSamplesSelected();",
        "        return;",
        "      }",
        "      setSelectedRuns(runs, 'page');",
        "    });",
        "    $removeSubset.on('click', function() {",
        "      selectedRuns = new Set();",
        "      applySelectionToCurrentPage();",
        "      sendSubsetRuns('reset', []);",
        "    });",
        "    $actions.append($selectedOnly, $pageOnly, $removeSubset);",
        "  }",
        "  function applySelectionToCurrentPage() {",
        "    syncingSelection = true;",
        "    table.rows({page: 'current'}).every(function() {",
        "      var run = runFromRowApi(this);",
        "      setRowSelected(this, selectedRuns.has(run));",
        "    });",
        "    syncingSelection = false;",
        "  }",
        "  $container.on('keyup change', 'div.dataTables_filter input', sendFilters);",
        "  $container.on('keyup change', 'thead tr input', sendFilters);",
        "  $container.on('click', 'tbody tr', function() {",
        "    rememberClickedRowSelection(this);",
        "  });",
        "  table.on('draw.dt', function() {",
        "    normalizeFilterLayout();",
        "    applySelectionToCurrentPage();",
        "    if (suppressNextFilterSend) {",
        "      suppressNextFilterSend = false;",
        "      return;",
        "    }",
        "    sendFilters();",
        "  });",
        "  installActionButtons();",
        "  normalizeFilterLayout();",
        "  sendFilters();",
        "  Shiny.addCustomMessageHandler('librariesApplyTableFilters', function(message) {",
        "    if (!message || message.table_id !== baseId) return;",
        "    var global = message.global || '';",
        "    var columns = message.columns || [];",
        "    globalSearchInput().val(global);",
        "    var $inputs = tableHeaderInputs();",
        "    $inputs.each(function(i){ $(this).val(columns[i] || ''); });",
        "    var changed = table.search() !== global;",
        "    if (changed) table.search(global);",
        "    for (var idx = 0; idx < table.columns().count(); idx++) {",
        "      var value = columns[idx] || '';",
        "      var column = table.column(idx);",
        "      if (column.search() === value) continue;",
        "      column.search(value);",
        "      changed = true;",
        "    }",
        "    if (changed) {",
        "      suppressNextFilterSend = true;",
        "      table.draw();",
        "    }",
        "  });",
        "  Shiny.addCustomMessageHandler('librariesApplyTableSelection', function(message) {",
        "    if (!message || message.table_id !== baseId) return;",
        "    selectedRuns = new Set(message.selected_runs || []);",
        "    applySelectionToCurrentPage();",
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
      keep <- keep & observatory_apply_dt_column_filter(df[[col_index]], pattern)
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

observatory_filter_range <- function(values, search_string) {
  if (!grepl("[.]{3}", search_string) || length(parts <- strsplit(search_string, "[.]{3}")[[1]]) > 2) {
    stop("The range of a numeric / date / time column must be of length 2")
  }
  if (length(parts) == 1) {
    parts <- c(parts, "")
  }
  parts <- gsub("^\\s+|\\s+$", "", parts)
  lower <- parts[1]
  upper <- parts[2]

  lower_num <- if (lower == "") NA_real_ else as.numeric(lower)
  upper_num <- if (upper == "") NA_real_ else as.numeric(upper)

  if (is.na(lower_num) && lower != "") {
    return(rep(FALSE, length(values)))
  }
  if (is.na(upper_num) && upper != "") {
    return(rep(FALSE, length(values)))
  }

  if (lower == "") {
    return(values <= upper_num)
  }
  if (upper == "") {
    return(values >= lower_num)
  }
  values >= lower_num & values <= upper_num
}

observatory_apply_dt_column_filter <- function(column, pattern) {
  if (is.null(pattern) || !nzchar(pattern)) {
    return(rep(TRUE, length(column)))
  }

  if (is.numeric(column)) {
    return(observatory_filter_range(column, pattern))
  }

  grepl(
    pattern,
    as.character(column),
    ignore.case = TRUE
  )
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
  displayed_libraries_df <- observatory$displayed_libraries_df %||%
    shiny::reactive(observatory_format_libraries_df(libraries_df()))
  filtered_displayed_libraries_df <- observatory$filtered_displayed_libraries_df %||%
    shiny::reactive(observatory_format_libraries_df(filtered_libraries_df()))

  data_table_filters <- shiny::reactiveVal(list())
  last_active_selection_id <- shiny::reactiveVal(NULL)
  restoring_active_plot_selection <- shiny::reactiveVal(FALSE)
  updating_active_subset <- shiny::reactiveVal(FALSE)

  displayed_column_count <- function() {
    ncol(displayed_libraries_df())
  }

  normalize_column_filters <- function(columns) {
    normalized <- rep("", displayed_column_count())
    if (is.null(columns)) return(normalized)
    columns <- as.character(columns)
    max_cols <- min(length(normalized), length(columns))
    if (max_cols > 0) {
      normalized[seq_len(max_cols)] <- columns[seq_len(max_cols)]
    }
    normalized
  }

  table_filters_for_selection <- function(selection_id) {
    filters <- data_table_filters()[[selection_id]]
    if (is.null(filters)) {
      filters <- list(global = "", columns = normalize_column_filters(NULL))
    } else {
      filters$global <- filters$global %||% ""
      filters$columns <- normalize_column_filters(filters$columns)
    }
    filters
  }

  table_filters_have_terms <- function(filters) {
    nzchar(filters$global %||% "") || any(nzchar(filters$columns))
  }

  with_active_plot_restore <- function(expr) {
    restoring_active_plot_selection(TRUE)
    session$onFlushed(
      function() restoring_active_plot_selection(FALSE),
      once = TRUE
    )
    force(expr)
  }

  with_active_subset_update <- function(expr) {
    updating_active_subset(TRUE)
    session$onFlushed(
      function() updating_active_subset(FALSE),
      once = TRUE
    )
    force(expr)
  }

  restore_active_plot_selection <- function(selection_id) {
    plot_sel <- selected_libraries$all_selections$plot_selections()[[selection_id]]
    plot_sel <- observatory_normalize_plot_selection(plot_sel, libraries_df()$Run)
    if (identical(current_plot_selection(), plot_sel)) {
      return(invisible(NULL))
    }

    with_active_plot_restore(current_plot_selection(plot_sel))
    invisible(NULL)
  }

  apply_data_table_filters <- function(selection_id) {
    filters <- table_filters_for_selection(selection_id)
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

  apply_data_table_selection <- function(selection_id, selected_runs = NULL) {
    if (missing(selected_runs)) {
      selected_runs <- selected_libraries$all_selections$data_table_selections()[[selection_id]]
    }
    if (is.null(selected_runs)) selected_runs <- character()
    if (observatory_selection_is_all_merged(selected_runs, libraries_df()$Run)) {
      selected_runs <- character()
    }
    session$sendCustomMessage(
      "librariesApplyTableSelection",
      list(
        table_id = session$ns("libraries_data_table"),
        selected_runs = as.character(selected_runs)
      )
    )
    invisible(NULL)
  }

  clear_active_table_filters <- function(selection_id) {
    filters <- data_table_filters()
    filters[[selection_id]] <- list(
      global = "",
      columns = normalize_column_filters(NULL)
    )
    data_table_filters(filters)
    apply_data_table_filters(selection_id)
    invisible(NULL)
  }

  commit_active_selection_runs <- function(runs) {
    selection_id <- selected_libraries$active_selection_id()
    if (is.null(selection_id) || selection_id == "") return(invisible(NULL))

    runs <- unique(as.character(runs))
    runs <- runs[!is.na(runs) & nzchar(runs)]
    if (length(runs) == 0) return(invisible(NULL))

    plot_selections <- selected_libraries$all_selections$plot_selections()
    plot_selections[[selection_id]] <- runs
    selected_libraries$all_selections$plot_selections(plot_selections)

    data_table_selections <- selected_libraries$all_selections$data_table_selections()
    data_table_selections[[selection_id]] <- runs
    selected_libraries$all_selections$data_table_selections(data_table_selections)

    data_table_selection_val(runs)
    selected_libraries$set_active_label(
      observatory_selection_subset_label(runs, libraries_df())
    )
    clear_active_table_filters(selection_id)
    with_active_subset_update(current_plot_selection(runs))
    session$sendCustomMessage("librariesActiveSelectionChanged", runs)
    apply_data_table_selection(selection_id)
    invisible(NULL)
  }

  reset_active_subset <- function() {
    selection_id <- selected_libraries$active_selection_id()
    if (is.null(selection_id) || selection_id == "") return(invisible(NULL))

    all_runs <- as.character(libraries_df()$Run)
    plot_selections <- selected_libraries$all_selections$plot_selections()
    plot_selections[[selection_id]] <- all_runs
    selected_libraries$all_selections$plot_selections(plot_selections)

    data_table_selections <- selected_libraries$all_selections$data_table_selections()
    data_table_selections[[selection_id]] <- all_runs
    selected_libraries$all_selections$data_table_selections(data_table_selections)

    data_table_selection_val(all_runs)
    selected_libraries$set_active_label("All merged")
    clear_active_table_filters(selection_id)
    with_active_subset_update(current_plot_selection(character()))
    session$sendCustomMessage("librariesActiveSelectionReset", "")
    apply_data_table_selection(selection_id, character())
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
    event <- input$libraries_data_table_subset_runs
    if (is.null(event)) return()
    if (identical(event$action, "reset")) {
      reset_active_subset()
    } else {
      commit_active_selection_runs(event$runs %||% character())
    }
  }) |> shiny::bindEvent(
    input$libraries_data_table_subset_runs,
    ignoreInit = TRUE
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
      filtered_libraries_df()
    ),
    {
      if (shiny::isolate(restoring_active_plot_selection())) {
        return()
      }
      if (shiny::isolate(updating_active_subset())) {
        return()
      }

      selection_value <- tryCatch(data_table_selection(), error = function(e) NULL)
      data_table_selection_val(selection_value)
    },
    ignoreInit = FALSE
  )

  shiny::observe({
    shiny::req(selected_libraries$active_selection_id())
    if (shiny::isolate(restoring_active_plot_selection())) {
      return()
    }
    if (shiny::isolate(updating_active_subset())) {
      return()
    }

    plot_sel <- current_plot_selection()
    if (!is.null(plot_sel) && length(plot_sel) != 0) {
      return()
    }

    selection_id <- selected_libraries$active_selection_id()
    filters <- data_table_filters()
    filters[[selection_id]] <- list(
      global = "",
      columns = normalize_column_filters(NULL)
    )
    data_table_filters(filters)
    selected_libraries$set_active_label("All merged")
    apply_data_table_filters(selection_id)
    apply_data_table_selection(selection_id, character())
  }) |> shiny::bindEvent(
    current_plot_selection(),
    ignoreInit = TRUE,
    ignoreNULL = FALSE
  )

  shiny::observe({
    shiny::req(selected_libraries$active_selection_id())
    if (shiny::isolate(restoring_active_plot_selection())) {
      return()
    }
    if (shiny::isolate(updating_active_subset())) {
      return()
    }

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
        columns = normalize_column_filters(input$libraries_data_table_manual_search_columns),
        global = input$libraries_data_table_manual_search
      )
      data_table_filters(filters)
    }
    last_active_selection_id(selection_id)
    restore_active_plot_selection(selection_id)
    apply_data_table_filters(selection_id)
    apply_data_table_selection(selection_id)
  }) |> shiny::bindEvent(selected_libraries$active_selection_id())

  shiny::observe({
    columns <- input$libraries_data_table_manual_search_columns
    global <- input$libraries_data_table_manual_search
    plot_sel <- current_plot_selection()
    if (!is.null(plot_sel) && length(plot_sel) > 0) {
      return()
    }
    if (is.null(columns) && is.null(global)) {
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
    input$libraries_data_table_manual_search
  )

  shiny::observe({
    shiny::req(selected_libraries$active_selection_id())
    selection_id <- selected_libraries$active_selection_id()
    filters <- data_table_filters()
    filters[[selection_id]] <- list(
      columns = normalize_column_filters(input$libraries_data_table_manual_search_columns),
      global = input$libraries_data_table_manual_search
    )
    data_table_filters(filters)
  }) |> shiny::bindEvent(
    input$libraries_data_table_manual_search,
    input$libraries_data_table_manual_search_columns
  )

  shiny::observe({
    shiny::req(selected_libraries$active_selection_id())
    selection_id <- selected_libraries$active_selection_id()
    filters <- table_filters_for_selection(selection_id)
    updated_df <- filtered_displayed_libraries_df()
    DT::replaceData(
      libraries_data_table_proxy(),
      updated_df,
      rownames = FALSE,
      resetPaging = TRUE,
      clearSelection = "none"
    )
    if (table_filters_have_terms(filters)) {
      apply_data_table_filters(selection_id)
    }
    apply_data_table_selection(selection_id)
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
  observatory_url_state = shiny::reactiveVal(NULL),
  templates = NULL
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

    displayed_libraries_df <- shiny::reactive({
      observatory_format_libraries_df(libraries_df())
    })

    output$libraries_umap_plot <- plotly::renderPlotly({
      observatory_selector_umap_plot_shiny(
        observatory_module = observatory_module(),
        color_by = input$color_by,
        session = session,
        templates = templates
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
      observatory_selector_data_table_shiny(
        libraries_df(),
        table_id = session$ns("libraries_data_table")
      )
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

    filtered_displayed_libraries_df <- shiny::reactive({
      displayed_libs <- displayed_libraries_df()
      plot_sel <- current_plot_selection()
      if (is.null(plot_sel) || length(plot_sel) == 0) {
        return(displayed_libs)
      }
      displayed_libs[Run %in% plot_sel]
    }) |> shiny::bindEvent(
      displayed_libraries_df(),
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
      displayed_column_count <- ncol(displayed_libraries_df())
      displayed_columns <- rep("", displayed_column_count)
      if (!is.null(manual_columns)) {
        manual_columns <- as.character(manual_columns)
        max_cols <- min(displayed_column_count, length(manual_columns))
        if (max_cols > 0) {
          displayed_columns[seq_len(max_cols)] <- manual_columns[seq_len(max_cols)]
        }
      }
      list(
        global = manual_global %||% "",
        columns = displayed_columns
      )
    })

    data_table_filtered_df <- shiny::reactive({
      filters <- current_table_filters()
      observatory_apply_dt_filters(
        filtered_displayed_libraries_df(),
        filters$global,
        filters$columns
      )
    })

    data_table_selection <- shiny::reactive({
      filters <- current_table_filters()
      plot_sel <- current_plot_selection()
      has_filters <- nzchar(filters$global) || any(nzchar(filters$columns))
      has_plot_subset <- !is.null(plot_sel) && length(plot_sel) > 0

      if (!has_filters && !has_plot_subset) return(libraries_df()$Run)
      if (!has_filters) return(filtered_libraries_df()$Run)
      data_table_filtered_df()$Run
    })

    observatory_selector_additional_controller(
      input, output, session,
      observatory = list(
        current_plot_selection = current_plot_selection,
        observatory_module = observatory_module,
        libraries_df = libraries_df,
        displayed_libraries_df = displayed_libraries_df,
        selected_libraries = selected_libraries,
        libraries_data_table_proxy = libraries_data_table_proxy,
        filtered_libraries_df = filtered_libraries_df,
        filtered_displayed_libraries_df = filtered_displayed_libraries_df,
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
