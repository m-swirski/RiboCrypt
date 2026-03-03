umap_ui <- function(id, all_exp_translons, label = "umap") {
  ns <- NS(id)
  tabPanel(
    title = "UMAP", icon = icon("rectangle-list"),
    h2("UMAP from count tables of all samples of all genes in organism"),
    # Include shinyjs so we can trigger hidden buttons
    fluidRow(
      column(2, umap_plot_type(ns)),
      column(2, experiment_input_select(all_exp_translons$name, ns)),
      column(2, umap_color_by_input_select(ns)),
      column(1, plot_button(ns("go")))
    ),
    tags$hr(),
    fluidRow(
      plotlyOutput(ns("c"), height = "700px") %>% shinycssloaders::withSpinner(color="#0dc5c1")
    )
  )
}


umap_server <- function(id, all_exp_meta, browser_options) {
  moduleServer(
    id,
    function(input, output, session) {
      plot_triggered <- reactiveVal(FALSE)
      md <- reactiveVal(NULL)
      default <- "all_samples-Homo_sapiens"
      selected_exp <- ifelse(default %in% all_exp_meta$name, default, "AUTO")
      experiment_update_select(NULL, all_exp_meta, all_exp_meta$name, selected_exp)
      # observeEvent(TRUE, {
      #   experiment_update_select(org, all_exp, experiments, rv$exp)
      # }, once = TRUE)

      # Trigger data loading when "Plot" is clicked.
      # Also reset the downloaded files vector.
      observeEvent(input$go, {
        md(load_data_umap(isolate(input$dff), isolate(input$umap_col)))
        plot_triggered(TRUE)
      })

      # Render DT Table ONLY if "Plot" was clicked
      output$c <- renderPlotly({
        req(plot_triggered())
        if (isolate(input$umap_plot_type) == "UMAP") {
          umap_plot(isolate(md()))
        } else umap_centroids_plot(isolate(md()))
        }) %>%
        bindCache(input$dff, input$umap_col, input$umap_plot_type) %>%
        bindEvent(input$go, ignoreInit = FALSE, ignoreNULL = TRUE)
    }
  )
}

umap_plot <- function(dt_umap, color.by = attr(dt_umap, "color.by")) {
  color.by.temp <- "color_column"
  names(color.by.temp) <- stringr::str_to_title(paste(color.by, collapse = " | "))
  color.by <- color.by.temp
  groups_n <- data.table::uniqueN(dt_umap$color_column)
  custom_palette <- grDevices::colorRampPalette(c(
    "#7A0403", # lava red
    "#FF6A00", # orange
    "#8B4513", # brown
    "#5E8C31", # olive green
    "#2E8B57", # sea green
    "#007FFF", # azure blue
    "#1F4FB2", # deep blue
    "#6A0DAD", # purple
    "#D65DB1", # pink
    "#F6B1B1"  # light red
  ))(max(groups_n, 3))

  dt_umap[, hover_text := paste0(
    "Bioproject ", BioProject, "<br />",
    "Run ID: ", sample, "<br />",
    "Author: ", author, "<br />",
    "Inhibitor: ", inhibitors, "<br />",
    "Tissue | CellLine: ", tissues_cell_lines
  )]

  plotly::plot_ly(
    data = dt_umap,
    x = ~`UMAP 1`,
    y = ~`UMAP 2`,
    type = "scatter",
    mode = "markers",
    color = ~color_column,
    colors = custom_palette,
    marker = list(size = 9, opacity = 0.75),
    text = ~hover_text,
    hovertemplate = "%{text}<extra></extra>"
  ) %>%
    plotly::layout(
      xaxis = list(title = "UMAP 1", zeroline = FALSE),
      yaxis = list(title = "UMAP 2", zeroline = FALSE),
      legend = list(title = list(text = names(color.by))),
      template = "plotly_white",
      dragmode = "select"
    ) |>
    plotly::config(
      displaylogo = FALSE,
      modeBarButtons = list(list("select2d", "lasso2d"))
    )
}

umap_centroids_plot <- function(dt_umap) {
  centroids <- dt_umap[, .(centroid_x = mean(`UMAP 1`), centroid_y = mean(`UMAP 2`)), by = color_column]
  coords <- as.matrix(centroids[, .(centroid_x, centroid_y)])
  rownames(coords) <- centroids$color_column
  dmat <- as.matrix(dist(coords))
  dt_dist <- as.data.table(as.table(dmat))
  setnames(dt_dist, c("Group1", "Group2", "Distance"))
  dt_dist <- dt_dist[Group1 != Group2]
  dt_dist[, Comparison := paste(Group1, "vs", Group2)]
  dt_dist <- dt_dist[as.integer(factor(Group1)) < as.integer(factor(Group2))]
  plot_ly(
    data = dt_dist,
    x = ~Group2,
    y = ~Group1,
    z = ~Distance,
    type = "heatmap",
    colorscale = "Viridis",
    hovertemplate = "Group1: %{y}<br>Group2: %{x}<br>Distance: %{z}<extra></extra>"
  ) %>%
    layout(
      title = "Upper Triangle: Pairwise Distance Between Group Centroids",
      xaxis = list(title = "Group 2", tickangle = 45),
      yaxis = list(title = "Group 1", autorange = "reversed")
    )
}
