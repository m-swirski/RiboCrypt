#' Compact, display-only controls for MegaBrowser.
#' @noRd
megabrowser_view_controls <- function(ns) {
  div(class = "mega-view-toolbar", div(class = "mega-view-controls",
      checkboxInput(ns("collapsed_clusters"), "Collapsed clusters", FALSE),
      selectizeInput(ns("visible_groups"), "Visible groups", choices = NULL, multiple = TRUE,
                     options = list(placeholder = "All groups")),
      selectizeInput(ns("enrichment_metadata"), "Enrichment metadata", choices = c("Grouping" = "grouping"))),
      div(class = "mega-layout-controls",
      tags$label(class = "mega-display-option", tags$input(type = "checkbox", checked = "checked",
        class = "mega-sidebar-control"), " Metadata sidebar"),
      tags$label(class = "mega-display-option", tags$input(id = ns("heatmap_only"), type = "checkbox",
        class = "mega-focus-control"), " Heatmap only"),
      tags$label(class = "mega-height-option", "Plot height",
        tags$input(type = "range", class = "mega-height-control", min = 400, max = 1200,
                   step = 100, value = 700, `aria-label` = "Plot height")),
      tags$button(id = ns("reset_view"), type = "button", class = "btn btn-outline-secondary mega-reset-view", title = "Reset zoom",
                  `aria-label` = "Reset zoom", icon("expand")),
      textOutput(ns("display_counts"), inline = TRUE)))
}

#' Distinguish focused libraries from displayed mean rows.
#' @noRd
megabrowser_display_counts <- function(original, display) {
  total <- ncol(original)
  visible <- if (isTRUE(attr(display$table, "collapsed_clusters"))) sum(display$meta$libraries) else ncol(display$table)
  count <- format(visible, big.mark = ",")
  if (visible != total) count <- paste(count, "of", format(total, big.mark = ","))
  sprintf("%s libraries | %s %s | %s positions", count,
          format(ncol(display$table), big.mark = ","), if (ncol(display$table) == 1L) "row" else "rows",
          format(nrow(original), big.mark = ","))
}

#' Scoped responsive MegaBrowser layout and client-only view preferences.
#' @noRd
megabrowser_layout_assets <- function() {
  tagList(tags$style(HTML(paste(
    ".mega-workspace {--mega-height:700px;}",
    ".mega-view-controls,.mega-layout-controls {display:flex;align-items:center;gap:20px;flex-wrap:wrap;padding:6px 0;}",
    ".mega-view-controls .form-group {margin:0;}",
    ".mega-view-controls .shiny-input-container {width:auto;}",
    ".mega-view-controls > .shiny-input-container:nth-child(2) {width:320px;}",
    ".mega-view-controls .checkbox {margin:0;}",
    ".mega-view-controls > .shiny-input-container:last-child {width:240px;margin-left:auto;}",
    ".mega-display-option,.mega-height-option {font-weight:normal;margin:0;}",
    ".mega-height-option {display:flex;align-items:center;gap:10px;}",
    ".mega-height-control {width:120px!important;}",
    ".mega-layout-controls .shiny-text-output {font-size:13px;color:#535c65;margin-left:auto;}",
    ".mega-focus-heatmap:not(.mega-static-view) .mega-plot-grid {grid-template-rows:0 100% 0;}",
    ".mega-focus-heatmap:not(.mega-static-view) .mega-summary-track,.mega-focus-heatmap:not(.mega-static-view) .mega-gene-track {display:none;}",
    ".mega-focus-heatmap:not(.mega-static-view) [id$='-myPlotlyPlot'] {height:var(--mega-height)!important;}",
    ".mega-plot-grid {display:grid;grid-template-columns:8% 92%;grid-template-rows:15% 75% 10%;height:var(--mega-height);width:100%;gap:0;}",
    ".mega-plot-sidebar {grid-column:1;grid-row:2;overflow:visible;}",
    ".mega-plot-main {grid-column:2;grid-row:1 / span 3;min-width:0;}",
    ".mega-workspace.mega-hide-sidebar .mega-plot-grid {grid-template-columns:0 100%;}",
    ".mega-workspace.mega-hide-sidebar .mega-plot-sidebar {display:none;}",
    ".mega-workspace .floating_settings_panel {width:min(680px,calc(100vw - 44px))!important;min-width:0!important;max-height:75vh;overflow:auto;}",
    ".mega-workspace .selectize-input {min-height:36px;height:auto;}",
    ".mega-workspace .selectize-input .item {overflow-wrap:anywhere;white-space:normal;}",
    "@media(max-width:767px){.mega-view-controls,.mega-layout-controls {gap:12px;}.mega-view-controls > .shiny-input-container:nth-child(2),.mega-view-controls > .shiny-input-container:last-child {margin-left:0;width:100%;}.mega-layout-controls .shiny-text-output {margin-left:0;}}"
  ))), tags$script(HTML(fetchJS("megabrowser_layout.js"))))
}
