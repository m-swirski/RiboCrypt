#' Lazy coding-region comparison workspace.
#' @noRd
megabrowser_translon_ui <- function(ns) {
  tabPanel("Translon enrichment", value = "Translon enrichment",
    div(class = "mega-layout-controls",
      selectizeInput(ns("translon_pair"), "Coverage ratio", choices = NULL),
      tags$span(icon("info-circle"), role = "img", `aria-label` = "Statistical interpretation",
        title = "Raw mean coverage per nucleotide; no pseudocount. Zero ratios are counted in statistics but omitted from log2 boxes. Rank tests are exploratory library-level comparisons: clusters use the same coverage, and libraries are not independent biological replicates."),
      downloadButton(ns("translon_ratios_csv"), "Library ratios CSV"),
      downloadButton(ns("translon_stats_csv"), "Statistics CSV"),
      downloadButton(ns("translon_regions_csv"), "Regions CSV")),
    tabsetPanel(tabPanel("Ratios", plotlyOutput(ns("translon_ratios_plot"))),
      tabPanel("Cluster statistics", DT::DTOutput(ns("translon_statistics"))),
      tabPanel("Regions",
        actionButton(ns("translon_add_region"), label = NULL, icon = icon("plus"),
          title = "Add transcript region", `aria-label` = "Add transcript region"),
        DT::DTOutput(ns("translon_regions")))))
}

#' Log ratios omit zeros/undefined values; the statistics table retains their counts.
#' @noRd
megabrowser_translon_ratio_plot <- function(ratios) {
  finite <- ratios[is.finite(Log2_ratio)]
  validate(need(nrow(finite) > 0L, "No positive, finite ratios for this pair. See Cluster statistics for zero and undefined counts."))
  plotly::plot_ly(finite, x = ~Group, y = ~Log2_ratio, type = "box", boxpoints = FALSE,
    text = ~Run, hoverinfo = "y", source = "mb_translons") %>%
    plotly::layout(xaxis = list(title = "Group", type = "category", categoryorder = "array", categoryarray = unique(ratios$Group)), yaxis = list(title = "log2 coverage-density ratio"),
      shapes = list(list(type = "line", xref = "paper", x0 = 0, x1 = 1, y0 = 0, y1 = 0,
                         line = list(dash = "dot", color = "#777777"))))
}

#' Tab gate is outside the cache so background observers cannot start the analysis.
#' @noRd
megabrowser_translon_outputs <- function(input, output, session, controller, table, grouped, workspace = NULL) {
  if (is.null(workspace)) workspace <- megabrowser_translon_workspace(controller, table)
  custom <- workspace$custom
  raw <- workspace$raw
  cached <- reactive(megabrowser_translon_analysis(controller(), table(), grouped(), custom(), raw())) %>%
    bindCache("translon-ratios-v2", controller()$table_hash, controller()$table_plot_hash, custom())
  analysis <- reactive({req(identical(input$mb_tabs, "Translon enrichment")); cached()})
  megabrowser_user_region_events(input, session, custom, raw)
  observeEvent(analysis(), {
    result <- analysis()
    labels <- setNames(result$regions$Labels, result$regions$Region)
    choices <- setNames(result$pairs$Pair, paste(labels[result$pairs$Numerator], labels[result$pairs$Denominator], sep = " / "))
    selected <- isolate(input$translon_pair)
    if (!length(selected) || !selected %in% choices) selected <- choices[1]
    updateSelectizeInput(session, "translon_pair", choices = choices, selected = selected, server = TRUE)
  })
  selected_ratios <- reactive({result <- analysis(); req(input$translon_pair); result$ratios[Pair == input$translon_pair]})
  selected_stats <- reactive({result <- analysis(); req(input$translon_pair); result$statistics[Pair == input$translon_pair]})
  output$translon_ratios_plot <- renderPlotly(megabrowser_translon_ratio_plot(selected_ratios()))
  output$translon_statistics <- DT::renderDT(DT::datatable(selected_stats(), rownames = FALSE,
    options = list(scrollX = TRUE)) %>% DT::formatRound(c("Q25", "Median", "Q75", "Ratio > 1 (%)", "Rank_biserial"), 3), server = TRUE)
  output$translon_regions <- DT::renderDT(analysis()$regions, rownames = FALSE, server = TRUE,
    options = list(scrollX = TRUE))
  output$translon_ratios_csv <- megabrowser_csv_download(reactive(analysis()$ratios), "translon-ratios")
  output$translon_stats_csv <- megabrowser_csv_download(reactive(analysis()$statistics), "translon-statistics")
  output$translon_regions_csv <- megabrowser_csv_download(reactive(analysis()$regions), "translon-regions")
  invisible(analysis)
}

#' Shared lazy raw coverage and custom regions for enrichment and display.
#' @noRd
megabrowser_translon_workspace <- function(controller, table) {
  custom <- reactiveVal(list())
  observeEvent(list(controller()$table_hash, controller()$table_plot_hash), custom(list()), priority = 100)
  raw <- reactive(as.matrix(load_collection(controller()$table_path, columns = colnames(table()$table))))
  regions <- reactive({
    base <- megabrowser_translon_regions(megabrowser_translon_panel(controller()), names(controller()$annotation))
    result <- c(base, custom())
    names(result) <- paste0("R", seq_along(result))
    result
  })
  list(custom = custom, raw = raw, regions = regions)
}

#' Modal errors remain visible without discarding the entered values.
#' @noRd
megabrowser_user_region_events <- function(input, session, custom, raw) {
  observeEvent(input$translon_add_region, {
    showModal(modalDialog(title = "Add transcript region",
      textInput(session$ns("translon_region_label"), "Label", paste0("UR", length(custom()) + 1L)),
      textInput(session$ns("translon_region_coordinates"), "Transcript coordinates", placeholder = "40:80"),
      uiOutput(session$ns("translon_region_error")),
      footer = tagList(modalButton("Cancel"), actionButton(session$ns("translon_region_save"), "Add region"))))
  })
  error <- reactiveVal(NULL)
  output <- session$output
  output$translon_region_error <- renderUI({req(error()); div(class = "text-danger", role = "alert", error())})
  observeEvent(input$translon_add_region, error(NULL))
  observeEvent(input$translon_region_save, {
    region <- tryCatch(megabrowser_user_region(input$translon_region_coordinates,
      input$translon_region_label, nrow(raw())), error = function(e) {error(conditionMessage(e)); NULL})
    if (is.null(region)) return()
    custom(append(custom(), list(region)))
    removeModal()
  })
}
