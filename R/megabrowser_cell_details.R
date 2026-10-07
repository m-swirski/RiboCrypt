#' Resolve heatmap coordinates through the saved row order, including group focus.
#' @noRd
megabrowser_cell_selection <- function(event, view, display, meta) {
  if (is.null(event) || !isTRUE(attr(display$table, "collapsed_translons"))) return(NULL)
  x <- suppressWarnings(as.numeric(event$x)); y <- suppressWarnings(as.numeric(event$y))
  if (length(x) != 1L || length(y) != 1L || !is.finite(x) || !is.finite(y)) return(NULL)
  if (x < 0.5 || x > nrow(view) + 0.5 || y < 0.5 || y > ncol(display$table) + 0.5) return(NULL)
  region <- floor(x + 0.5); row <- floor(y + 0.5)
  if (region > nrow(view) || row > ncol(display$table)) return(NULL)
  orders <- unlist(attr(display$table, "row_order_list"), use.names = FALSE)
  column <- colnames(display$table)[orders[row]]
  list(region = rownames(view)[region], column = column,
       runs = colnames(view)[megabrowser_cell_memberships(view, display, meta)[[column]]])
}

#' Full analysis libraries are reused for ratio-by-metadata comparisons.
#' @noRd
megabrowser_cell_libraries <- function(selection, view, meta, metadata) {
  runs <- colnames(view)
  result <- copy(meta[match(runs, Run)])
  for (field in setdiff(names(metadata), names(result))) set(result, j = field, value = metadata[[field]][match(runs, metadata$Run)])
  region <- attr(view, "translon_regions")[Region == selection$region]
  result[, `:=`(Run = runs, Selected = runs %in% selection$runs,
    Density = as.numeric(view[selection$region, ]), Coverage_sum = as.numeric(view[selection$region, ]) * region$Length_nt)]
  result
}

#' Zero numerators stay valid; zero/undefined reference densities stay undefined.
#' @noRd
megabrowser_cell_ratios <- function(libraries, view, reference, minimum = 0) {
  req(length(reference) == 1L, reference %in% rownames(view))
  validate(need(length(minimum) == 1L && is.finite(minimum) && minimum >= 0,
    "Minimum coverage sum must be a finite non-negative number."))
  result <- copy(libraries)
  denominator <- as.numeric(view[reference, match(result$Run, colnames(view))])
  reference_sum <- denominator * attr(view, "translon_regions")[Region == reference, Length_nt]
  result[, `:=`(Reference_density = denominator,
    Reference_coverage_sum = reference_sum,
    Supported = is.finite(Coverage_sum) & Coverage_sum > 0 & Coverage_sum >= minimum &
      is.finite(reference_sum) & reference_sum > 0 & reference_sum >= minimum,
    Ratio = ifelse(is.finite(denominator) & denominator > 0 & is.finite(Density), Density / denominator, NA_real_))]
  result[, Log2FC := log2(Ratio)]
  result
}

#' Cheap term-versus-rest summaries are exploratory, not study-adjusted effects.
#' @noRd
megabrowser_ratio_metadata <- function(libraries, field) {
  req(field %in% names(libraries))
  data <- copy(libraries)
  values <- data[[field]]
  data[, Term := if (is.numeric(values)) as.character(megabrowser_numeric_bins(values)) else
    fifelse(is.na(values) | as.character(values) == "", "(Missing)", as.character(values))]
  study <- intersect(c("study", "bioproject", "project", "study_accession", "project_accession", "study_id"), tolower(names(data)))
  data[, Study := if (length(study)) as.character(data[[match(study[1], tolower(names(data)))]]) else NA_character_]
  result <- data[, {
    ratios <- Ratio[is.finite(Ratio)]
    rest <- data[Term != .BY$Term & is.finite(Ratio), Ratio]
    q <- if (length(ratios)) unname(quantile(ratios, c(0.25, 0.5, 0.75))) else rep(NA_real_, 3)
    studies <- Study[!is.na(Study) & nzchar(Study)]
    c(list(Libraries = .N, Selected = sum(Selected), Supported = if ("Supported" %in% names(data)) sum(Supported) else NA_integer_,
      Valid = length(ratios), Undefined = sum(!is.finite(Ratio)),
      Zero = sum(ratios == 0), Q25 = log2(q[1]), Median = log2(q[2]), Q75 = log2(q[3]),
      Studies = if (length(study)) uniqueN(studies) else NA_integer_,
      `Largest study (%)` = if (length(studies)) 100 * max(table(studies)) / .N else NA_real_),
      megabrowser_translon_test(ratios, rest))
  }, by = Term]
  result[, BH := p.adjust(P, "BH")]
  result <- result[order(-Median, Term, na.last = TRUE)]
  setnames(result, c("Q25", "Median", "Q75"), c("Q25_log2FC", "Median_log2FC", "Q75_log2FC"))
  result
}

#' Plot finite log2 ratios; zero ratios remain explicit in counts and exports.
#' @noRd
megabrowser_cell_distribution <- function(libraries, metric = "Density") {
  data <- copy(libraries)
  data[, Value := get(metric)]
  if (metric == "Ratio") data[, Value := log2(Value)]
  data <- data[is.finite(Value)]
  validate(need(nrow(data) > 0, "No finite values; log2 scores require positive region and reference coverage."))
  data[, Set := ifelse(Selected, "Selected", "Other libraries")]
  plotly::plot_ly(data, x = ~Set, y = ~Value, type = "box", boxpoints = FALSE) %>%
    plotly::layout(xaxis = list(title = NULL), yaxis = list(title = if (metric == "Ratio") "log2(region / reference density)" else "Mean coverage per nucleotide"))
}

#' Cell inspector is opened only by an explicit heatmap click.
#' @noRd
megabrowser_cell_modal <- function(ns, selection, regions, fields, field) {
  region <- regions[Region == selection$region]
  references <- regions[Type %in% c("clean_cds", "CDS")]
  references <- references[order(Type != "clean_cds")]
  modalDialog(title = paste(region$Display_label, "-", selection$column), size = "l",
    tags$p(style = "overflow-wrap:anywhere", paste(region$Region, region$Labels, "|", region$Intervals, "|", region$Length_nt, "nt")),
    uiOutput(ns("cell_summary")),
    div(class = "mega-layout-controls",
      selectInput(ns("cell_metric"), "Distribution", c("Raw coverage" = "Density", "Region / reference (log2 FC)" = "Ratio")),
      selectInput(ns("cell_reference"), "Reference region", setNames(references$Region, paste(references$Display_label, references$Region))),
      selectizeInput(ns("cell_metadata"), "Metadata field", fields, selected = field),
      tags$span(icon("info-circle"), title = "Scores are log2(region/reference density), without a pseudocount. Zero ratios have score -Inf and are counted but omitted from the finite distribution. Metadata quartiles are log2-transformed raw-ratio quartiles, including zeros. Support requires adequate numerator and reference coverage but does not filter rank tests. Coverage sums are not independent read counts. Tests use all analysis libraries, retain zeros and exclude undefined denominators, with BH adjustment within the chosen field; not study-adjusted or replicate-aware evidence.")),
    tabsetPanel(id = ns("cell_tabs"), tabPanel("Distribution", plotlyOutput(ns("cell_distribution"), height = "260px")),
      tabPanel("Libraries", DT::DTOutput(ns("cell_libraries"))),
      tabPanel("Ratio by metadata", DT::DTOutput(ns("cell_ratio_metadata")))),
    footer = tagList(downloadButton(ns("cell_download"), "Libraries CSV"), modalButton("Close")))
}

#' Dynamic controls must use the same Bootstrap theme as the app shell.
#' @noRd
megabrowser_show_cell_modal <- function(session, cell, regions, fields, field) {
  current <- session$getCurrentTheme()
  if (is.null(current)) current <- rc_theme()
  shiny::withReactiveDomain(NULL, {
    theme <- shiny::getShinyOption("bootstrapTheme")
    shiny::shinyOptions(bootstrapTheme = current)
    tryCatch(showModal(megabrowser_cell_modal(session$ns, cell, regions, fields, field), session = session),
      finally = shiny::shinyOptions(bootstrapTheme = theme))
  })
}

#' Reset selection before a new display can reinterpret old plot coordinates.
#' @noRd
megabrowser_cell_outputs <- function(input, output, session, view, display, grouped, metadata) {
  selected <- reactiveVal(NULL)
  observeEvent(display(), {if (!is.null(selected())) {selected(NULL); removeModal()}}, ignoreInit = TRUE, priority = 100)
  observeEvent(get_plotly_session_event(session, "plotly_click", "mb_mid"), {
    event <- get_plotly_session_event(session, "plotly_click", "mb_mid")
    cell <- megabrowser_cell_selection(event, view(), display(), grouped()$meta)
    if (is.null(cell)) return()
    selected(cell)
    fields <- unique(c(setdiff(names(grouped()$meta), c("Run", "index", "order", "cluster")), setdiff(names(metadata), "Run")))
    field <- if (isTruthy(input$enrichment_metadata) && input$enrichment_metadata %in% fields) input$enrichment_metadata else fields[1]
    megabrowser_show_cell_modal(session, cell, attr(view(), "translon_regions"), fields, field)
  })
  libraries <- reactive({req(selected()); megabrowser_cell_libraries(selected(), view(), grouped()$meta, metadata)})
  ratios <- reactive(megabrowser_cell_ratios(libraries(), view(), input$cell_reference, input$region_min_coverage))
  output$cell_summary <- renderUI({
    data <- libraries()[Selected == TRUE]
    min_cov <- input$region_min_coverage
    supported <- sum(is.finite(data$Coverage_sum) & data$Coverage_sum > 0 & data$Coverage_sum >= min_cov)
    score_counts <- if (identical(input$cell_metric, "Ratio")) {
      scores <- ratios()[Selected == TRUE]
      paste("|", sum(scores$Ratio == 0, na.rm = TRUE), "zero ratios (-Inf);",
        sum(!is.finite(scores$Ratio)), "undefined reference ratios")
    } else ""
    tags$p(paste(nrow(data), "libraries |", supported, "meet coverage sum >=", min_cov,
      "| Mean density:", signif(mean(data$Density), 4), score_counts))
  })
  output$cell_distribution <- renderPlotly({
    req(input$cell_metric)
    megabrowser_cell_distribution(if (input$cell_metric == "Ratio") ratios() else libraries(), input$cell_metric)
  })
  output$cell_libraries <- DT::renderDT({
    data <- if (isTruthy(input$cell_reference)) ratios() else libraries()
    fields <- intersect(c("Run", "Density", "Coverage_sum", "Reference_density", "Reference_coverage_sum", "Log2FC", "Ratio", "Supported", "grouping", "TISSUE", "CELL_LINE", "CONDITION", "BioProject", "study"), names(data))
    DT::datatable(data[Selected == TRUE, ..fields], rownames = FALSE, options = list(scrollX = TRUE, pageLength = 10))
  }, server = TRUE)
  output$cell_ratio_metadata <- DT::renderDT({
    req(identical(input$cell_tabs, "Ratio by metadata"))
    DT::datatable(megabrowser_ratio_metadata(ratios(), input$cell_metadata), rownames = FALSE,
      options = list(scrollX = TRUE, pageLength = 10)) %>% DT::formatRound(c("Q25_log2FC", "Median_log2FC", "Q75_log2FC", "Largest study (%)", "Rank_biserial"), 3)
  }, server = TRUE)
  output$cell_download <- megabrowser_csv_download(reactive({data <- if (isTruthy(input$cell_reference)) ratios() else libraries(); data[Selected == TRUE]}), "region-cell-libraries")
  invisible(selected)
}
