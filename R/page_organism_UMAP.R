umap_ui <- function(id, all_exp_translons, label = "umap") {
  ns <- NS(id)
  tabPanel(
    title = "UMAP", icon = icon("rectangle-list"),
    h2("UMAP from count tables of all samples of all genes in organism"),
    # Include shinyjs so we can trigger hidden buttons
    sidebarLayout(
      sidebarPanel(
        experiment_input_select(all_exp_translons$name, ns),
        actionButton(ns("go"), "Search", icon = icon("magnifying-glass"))
      ),
      mainPanel(
        plotlyOutput(ns("c"), height = "700px") %>% shinycssloaders::withSpinner(color="#0dc5c1")
      )
    )
  )
}


umap_server <- function(id) {
  moduleServer(
    id,
    function(input, output, session) {
      plot_triggered <- reactiveVal(FALSE)
      md <- reactiveVal(NULL)
      # Trigger data loading when "Plot" is clicked.
      # Also reset the downloaded files vector.
      observeEvent(input$go, {
        md(load_data_umap(isolate(input$dff)))
        plot_triggered(TRUE)
      })

      # Render DT Table ONLY if "Plot" was clicked
      output$c <- renderPlotly({req(plot_triggered()); umap_plot(isolate(md()$dt_umap))}) %>%
        bindCache(input$dff) %>%
        bindEvent(input$go, ignoreInit = FALSE, ignoreNULL = TRUE)

      check_url_for_basic_parameters()
    }
  )
}

umap_plot <- function(dt_umap) {
  gg <- ggplot(dt_umap, aes(x = `UMAP 1`, y = `UMAP 2`, color = tissues_cell_lines)) + geom_point() + cowplot::theme_cowplot() + scale_fill_viridis_b()
  text_aes <- aes(text = paste0(
    "Bioproject ", dt_umap$BioProject, "\n",
    "Run ID: ",dt_umap$sample, "\n",
    "Author: ", dt_umap$author, "\n",
    "Inhibitor: ", dt_umap$inhibitors, "\n",
    "Tissue | CellLine: ", dt_umap$tissues_cell_lines
  ))
  plotly::ggplotly(gg + text_aes, tooltip = "text")
}

