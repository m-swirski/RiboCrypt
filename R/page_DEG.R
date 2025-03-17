DEG_ui <- function(id, all_exp, browser_options, label = "DEG") {
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
                   experiment_input_select(experiments, ns, browser_options),
                   diff_method_input_select(ns),
                   condition_input_select(ns),
                   library_input_select(ns, label = "Libraries (Group1)", id = "library1"),
                   library_input_select(ns, label = "Libraries (Group2)", id = "library2"),
                   helper_button_redirect_call()
          ),
          tabPanel("Settings",
                   checkboxInput(ns("draw_unnreg"), label = "Draw unregulated", value = FALSE),
                   checkboxInput(ns("other_tx"), label = "Full annotation (all isoforms)", value = FALSE),
                   sliderInput(ns("pval"), "P-value", min = 0, max = 1,
                               value = 0.05, step = 0.01),
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


DEG_server <- function(id, all_experiments, env, df, experiments, libs,
                       org, rv) {
  moduleServer(
    id,
    function(input, output, session, all_exp = all_experiments) {

      # Update main side panels
      uses_gene <- FALSE
      study_and_gene_observers(input, output, session)
      cond <- reactive(
        if (nrow(df()) > 1) {
          design <- design(df(), batch.correction.design = TRUE, multi.factor = FALSE)
          target.contrast <- design[1]
          df()[, target.contrast]
        } else "")
      observeEvent(cond(), condition_update_select(cond))

      observeEvent(cond(), {
        cond_1 <- cond()[1]
        library_update_select(libs, libs()[cond() == cond_1], "library1")
        }, ignoreNULL = TRUE, ignoreInit = TRUE)

      observeEvent(cond(), {
        cond_2 <- unique(cond())
        req(length(cond_2) > 1)
        cond_2 <- cond_2[2]
        library_update_select(libs, libs()[cond() == cond_2], "library2")
      }, ignoreNULL = TRUE, ignoreInit = TRUE)

      # Main plot, this code is only run if 'plot' is pressed
      controls <- eventReactive(input$go,
                                click_plot_DEG_main_controller(input, df, libs))
      model <- reactive(DE_model_from_ctrl(controls)) %>%
        bindCache(controls()$hash_string_pre)
      #
      analysis_dt <- reactive(DE_model_results(model(), controls)) %>%
        bindCache(controls()$hash_string_full)
      output$c <- renderPlotly(DEG_plot(analysis_dt(),
                draw_non_regulated = controls()$draw_unregulated)) %>%
        bindCache(controls()$hash_string_full)
      return(rv)
    }
  )
}
