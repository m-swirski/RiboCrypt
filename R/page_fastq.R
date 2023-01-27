fastq_ui <- function(id, label = "fastq", all_exp) {
  ns <- NS(id)
  genomes <- unique(all_exp$organism)
  experiments <- all_exp$name
  tabPanel(
    title = "fastq", icon = icon("layer-group"),
    sidebarLayout(
      sidebarPanel(
        tabPanel("fastq",
                 organism_input_select(c("ALL", genomes), ns),
                 experiment_input_select(experiments, ns),
                 library_input_select(ns, FALSE)),
        actionButton(ns("go"), "find", icon = icon("magnifying-glass")),
      ),
      mainPanel(
        htmlOutput(ns('fastq')) %>% shinycssloaders::withSpinner(color="#0dc5c1")
      )
    )
  )
}

fastq_server <- function(id, all_experiments, relative_dir_to_bam = "../trim",
                         env = environment()) {
  moduleServer(
    id,
    function(input, output, session, all_exp = all_experiments,
             relative_dir = relative_dir_to_bam) {
      # Organism / study objects
      org_and_study_changed_checker(input, output, session)
      # Update main side panels
      observeEvent(org(), experiment_update_select(org, all_exp, experiments))
      observeEvent(libs(), library_update_select(libs))

      # Main plot controller, this code is only run if 'plot' is pressed
      page <- eventReactive(input$go,
                 get_fastq_page(input, libs, df, relative_dir))
      # Main plot, this code is only run if 'plot' is pressed
      output$fastq <- renderUI({page()})
    }
  )
}
