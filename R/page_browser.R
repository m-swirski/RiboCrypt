browser_ui = function(id, label = "Browser", all_exp) {
  ns <- NS(id)
  genomes <- unique(all_exp$organism)
  experiments <- all_exp$name
  tabPanel(
    title = "browser", icon = icon("chart-line"),
    sidebarLayout(
      sidebarPanel(
        tabsetPanel(
          tabPanel("Browser",
                   tabPanel("a", tags$head(tags$script(HTML('
                        var fakeClick = function(tabName, anchorName) {
                          var dropdownList = document.getElementsByTagName("a");
                          for (var i = 0; i < dropdownList.length; i++) {
                            var link = dropdownList[i];
                            if(link.getAttribute("data-value") == tabName) {
                              link.click();
                              document.getElementById(anchorName).scrollIntoView({
                                behavior: "smooth"
                                });
                            };
                          }
                        };
                                                 ')))),
                   organism_input_select(c("ALL", genomes), ns),
                   experiment_input_select(experiments, ns),
                   gene_input_select(ns),
                   tx_input_select(ns),
                   library_input_select(ns),
                   frame_type_select(ns),
                   sliderInput(ns("kmer"), "K-mer length", min = 1, max = 20, value = 1)
          ),
          tabPanel("Settings",
                   numericInput(ns("extendLeaders"), "5' extension", 0),
                   numericInput(ns("extendTrailers"), "3' extension", 0),
                   checkboxInput(ns("viewMode"), label = "Genomic View", value = FALSE),
                   checkboxInput(ns("useCustomRegions"), label = "Protein structures", value = FALSE),
                   checkboxInput(ns("other_tx"), label = "Full annotation", value = FALSE),
                   checkboxInput(ns("summary_track"), label = "Summary top track", value = FALSE),
                   frame_type_select(ns, "summary_track_type", "Select summary display type"),
                   export_format_of_plot(ns)
          ),
        ),
         actionButton(ns("go"), "Plot", icon = icon("rocket")),
      ),
      mainPanel(
        plotlyOutput(outputId = ns("c")) %>% shinycssloaders::withSpinner(color="#0dc5c1"),
        uiOutput(ns("variableUi"))
      )
    )
  )
}

browser_server <- function(id, all_experiments, env) {
  moduleServer(
    id,
    function(input, output, session, all_exp = all_experiments) {
      # Static values
      genomes <- unique(all_exp$organism)
      experiments <- all_exp$name

      # Set reactive values
      org <- reactive(input$genome)
      rv <- reactiveValues(lstval="",curval="") # Store current and last genome
      rv_changed <- reactiveVal(NULL) # Did genome change?
      df <- reactive(get_exp(input$dff, experiments, env))
      observeEvent(df(), update_rv(rv, df), priority = 2)
      observe(update_rv_changed(rv, rv_changed), priority = 1) %>%
        bindEvent(rv$curval)
      tx <- reactive({loadRegion(df())}) %>%
        bindEvent(rv_changed(), ignoreNULL = TRUE)
      cds <- reactive(loadRegion(df(), part = "cds")) %>%
        bindEvent(rv_changed(), ignoreNULL = TRUE)
      gene_name_list <- reactive(get_gene_name_categories(df())) %>%
        bindCache(rv$curval) %>%
        bindEvent(rv_changed(), ignoreNULL = TRUE)
      libs <- reactive(bamVarName(df()))
      # Update main side panels
      observeEvent(org(), experiment_update_select(org, all_exp, experiments))
      observeEvent(gene_name_list(), gene_update_select(gene_name_list))
      observeEvent(input$gene, tx_update_select(isolate(input$gene),
                      gene_name_list), ignoreNULL = TRUE, ignoreInit = TRUE)
      observeEvent(libs(), library_update_select(libs))

      # Main plot controller, this code is only run if 'plot' is pressed
      mainPlotControls <- eventReactive(input$go,
        click_plot_browser_main_controller(input, tx, cds, libs, df))
      # Main plot, this code is only run if 'plot' is pressed
      output$c <- renderPlotly(click_plot_browser(mainPlotControls, session))
      # source(system.file("R", "Module_Protein_structures.R", package = "RiboCrypt"),  local = TRUE)
      module_protein(input, output, session)
    }
  )
}
