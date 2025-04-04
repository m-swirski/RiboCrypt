fastq_ui <- function(id, all_exp, browser_options, libs, label = "fastq") {
  ns <- NS(id)
  genomes <- unique(all_exp$organism)
  experiments <- all_exp$name
  tabPanel(
    title = "FastQ report", icon = icon("layer-group"),
    sidebarLayout(
      sidebarPanel(
        tabPanel("fastq",
                 organism_input_select(c("ALL", genomes), ns),
                 experiment_input_select(experiments, ns, browser_options),
                 library_input_select(ns, FALSE, libs)),
        plot_button(ns("go"), "View report", icon_type = icon("magnifying-glass")),
        width = 3),
      mainPanel(
        htmlOutput(ns('fastq')) %>% shinycssloaders::withSpinner(color="#0dc5c1")
      )
    )
  )
}

fastq_server <- function(id, all_experiments, df, experiments, libs, org, rv,
                         relative_dir_to_bam = "../trim",
                         env = environment()) {
  moduleServer(
    id,
    function(input, output, session, all_exp = all_experiments,
             relative_dir = relative_dir_to_bam) {
      uses_gene <- FALSE
      study_and_gene_observers(input, output, session)

      page <- reactiveVal(NULL)  # Store the path separately

      observeEvent(input$go, {
        path <- get_fastq_page(input, libs, df, relative_dir)
        if (!is.null(path)) {
          page(NULL)  # Reset first
          invalidateLater(700)  # Delay reactivity by 0.5s
          page(path)  # Update the path
        }
      })

      output$fastq <- renderUI({
        req(page())  # Ensure path exists
        tags$iframe(seamless = "seamless", src = page(), width = 1000, height = 900)
      })
      return(rv)
    }
  )
}
