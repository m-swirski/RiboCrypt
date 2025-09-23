umap_ui <- function(id, all_exp_translons, gene_names_init, browser_options, label = "umap") {
  ns <- NS(id)
  tabPanel(
    title = "UMAP", icon = icon("rectangle-list"),
    h2("UMAP from count tables of all samples of all genes in organism"),
    # Include shinyjs so we can trigger hidden buttons
    tags$hr(),
    tabsetPanel(
      tabPanel("UMAP plot",
               fluidRow(
                 column(2, umap_plot_type(ns)),
                 column(2, experiment_input_select(all_exp_translons$name, ns)),
                 column(2, umap_color_by_input_select(ns)),
                 column(1, plot_button(ns("go")))
               ),
               fluidRow(
                 plotlyOutput(ns("c"), height = "700px") %>% shinycssloaders::withSpinner(color="#0dc5c1")
               ))
    ),
    sampleSelectionsUi(ns("sampleSelection"))
  )
}


umap_server <- function(id, metadata, all_exp_meta, tx, cds, libs, df, browser_options) {
  moduleServer(
    id,
    function(input, output, session) {
      ns <- session$ns
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
      output$c <- {
        renderPlotly({
          generated_plot <- {
            req(plot_triggered())
            if (isolate(input$umap_plot_type) == "UMAP") {
              umap_plot(isolate(md()$dt_umap))
            } else umap_centroids_plot(isolate(md()$dt_umap))
          }
          onRender(generated_plot, fetchJS("umap_plot_extension.js"), ns("selectedPoints"))
        }) %>%
        bindCache(input$dff, input$umap_col, input$umap_plot_type) %>%
        bindEvent(input$go, ignoreInit = FALSE, ignoreNULL = TRUE)
      }
      
      
      selectedSamples <- sampleSelectionsServer("sampleSelection", metadata, reactive(input$selectedPoints))

      check_url_for_basic_parameters()
    }
  )
}

umap_plot <- function(dt_umap, color.by = attr(dt_umap, "color.by")) {
  color.by.temp <- "color_column"
  names(color.by.temp) <- stringr::str_to_title(paste(color.by, collapse = " | "))
  color.by <- color.by.temp

  gg <- ggplot(dt_umap, aes(x = `UMAP 1`, y = `UMAP 2`, color = color_column)) +
    geom_point() + cowplot::theme_cowplot() + scale_fill_viridis_b() +
    labs(color = names(color.by))
  text_aes <- aes(text = paste0(
    "Bioproject ", dt_umap$BioProject, "\n",
    "Run ID: ",dt_umap$sample, "\n",
    "Author: ", dt_umap$author, "\n",
    "Inhibitor: ", dt_umap$inhibitors, "\n",
    "Tissue | CellLine: ", dt_umap$tissues_cell_lines
  ))
  plotly::ggplotly(gg + text_aes, tooltip = "text")
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
