DEG_ui = function(id, label = "DEG", all_exp) {
  ns <- NS(id)
  genomes <- unique(all_exp$organism)
  experiments <- all_exp$name
  tabPanel(
    title = "Differential expression", icon = icon("layer-group"),
    sidebarLayout(
      jqui_resizable(jqui_draggable(sidebarPanel(
        tabsetPanel(
          tabPanel("Differential expression",
                   organism_input_select(c("ALL", genomes), ns),
                   experiment_input_select(experiments, ns),
                   library_input_select(ns),
                   condition_input_select(ns),
                   helper_button_redirect_call()
          ),
          tabPanel("Settings",
                  
                   checkboxInput(ns("other_tx"), label = "Full annotation", value = FALSE),
                   export_format_of_plot(ns)
          ),
        ),
        actionButton(ns("go"), "Plot", icon = icon("rocket")),
      ))),
      mainPanel(
        jqui_resizable(plotlyOutput(outputId = ns("c"), height = "500px")) %>% shinycssloaders::withSpinner(color="#0dc5c1"),
        uiOutput(ns("variableUi")
      ))
    )
  )
}


DEG_server <- function(id, all_experiments, env) {
  moduleServer(
    id,
    function(input, output, session, all_exp = all_experiments) {
      # Organism / study objects
      org_and_study_changed_checker(input, output, session)
      # Gene objects
      
      valid_genes_subset <- reactive(filterTranscripts(df(), stopOnEmpty = FALSE,
                                                       minFiveUTR = 0, minThreeUTR = 0)) %>%
        bindEvent(rv_changed(), ignoreNULL = TRUE)
      length_table <- reactive(optimizedTranscriptLengths(df(), TRUE, TRUE, FALSE)) %>%
        bindEvent(rv_changed(), ignoreNULL = TRUE)
      tx <- reactive({loadRegion(df(), part = "mrna", names.keep = valid_genes_subset())}) %>%
        bindEvent(rv_changed(), ignoreNULL = TRUE)
      cds <- reactive(loadRegion(df(), part = "cds", names.keep = valid_genes_subset())) %>%
        bindEvent(rv_changed(), ignoreNULL = TRUE)
      gene_name_list <- reactive(get_gene_name_categories(df())) %>%
        bindEvent(rv_changed(), ignoreNULL = TRUE)
      libs <- reactive(bamVarName(df()))
      # Update main side panels
      all_is_gene <- TRUE
      study_and_gene_observers(input, output, session)
      cond <- reactive(df()$condition)
      observeEvent(cond(), condition_update_select(cond))
      
      
      # Main plot, this code is only run if 'plot' is pressed
      
      mainPlotControls <- eventReactive(input$go, 
                                   
                                        click_plot_DEG_main_controller(input, df))
      analysis_dt <- reactive(DEG.analysis(mainPlotControls()$dff,  output.dir = NULL))
      output$c <- renderPlotly(DEG_plot(analysis_dt(), draw_non_regulated = TRUE))
      # output$c <- renderPlotly(ggplotly(ggplot(aes(x= 1:length(analysis_dt()), y = 1:length(analysis_dt()))) + geom_point()))

    }
  )
}