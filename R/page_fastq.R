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

fastq_server <- function(id, all_experiments, relative_dir_to_bam = "../trim") {
  moduleServer(
    id,
    function(input, output, session, all_exp = all_experiments,
             relative_dir = relative_dir_to_bam) {
      # Static values
      genomes <- unique(all_exp$organism)
      experiments <- all_exp$name

      # Set reactive values
      org <- reactive(input$genome)
      rv <- reactiveValues(lstval="",curval="") # Store current and last genome
      rv_changed <- reactiveVal(NULL) # Did genome change?
      df <- reactive(get_exp(input$dff, experiments, environment()))
      observeEvent(df(), update_rv(rv, df), priority = 2)
      observe(update_rv_changed(rv, rv_changed), priority = 1) %>%
        bindEvent(rv$curval)
      libs <- reactive(bamVarName(df()))
      # Update main side panels
      observeEvent(org(), experiment_update_select(org, all_exp, experiments))
      observeEvent(libs(), library_update_select(libs))

      # Main plot controller, this code is only run if 'plot' is pressed
      page <- eventReactive(input$go, {
          print("In fastq page")
          dff <- observed_exp_subset(isolate(input$library), libs, df)
          trim_dir <- file.path(libFolder(dff), relative_dir)
          if (!dir.exists(trim_dir)) {
            warning("No valid trim directory")
            return(NULL)
          }
          candidates <- list.files(trim_dir, full.names = TRUE, pattern = "html")
          candidates_base <- gsub(".html$", "", basename(candidates))
          proper_names <- gsub("_Aligned.*", "", ORFik:::remove.file_ext(dff$filepath,basename = TRUE))
          path <- grep(pattern = proper_names, candidates, value = TRUE)
          if (length(path) != 1) {
            hits <- lapply(candidates_base, function(x) grep(x, proper_names))
            path <- candidates[hits]
            if (length(path) != 1) {
              warning("No valid html file found in folder!")
              return(NULL)
            }
          }
          print(path)
          addResourcePath("tmpuser", dirname(path))
          path <- file.path("tmpuser", basename(path))
          page <- tags$iframe(seamless="seamless", src= path, width=1000, height=900)
        })
      # Main plot, this code is only run if 'plot' is pressed
      output$fastq <- renderUI({page()})
    }
  )
}
