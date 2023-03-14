fastq_ui <- function(id, all_exp, browser_options, libs, label = "fastq") {
  ns <- NS(id)
  genomes <- unique(all_exp$organism)
  experiments <- all_exp$name
  tabPanel(
    title = "fastq", icon = icon("layer-group"),
    sidebarLayout(
      sidebarPanel(
        tabPanel("fastq",
                 organism_input_select(c("ALL", genomes), ns),
                 experiment_input_select(experiments, ns, browser_options),
                 library_input_select(ns, FALSE, libs)),
        actionButton(ns("go"), "find", icon = icon("magnifying-glass")),
      ),
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
      # Main plot controller, this code is only run if 'plot' is pressed
      page <- eventReactive(input$go,
                 get_fastq_page(input, libs, df, relative_dir))
      # Main plot, this code is only run if 'plot' is pressed
      output$fastq <- renderUI({page()})
      return(rv)
    }
  )
}
