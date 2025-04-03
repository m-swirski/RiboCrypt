codon_ui <- function(id, all_exp, browser_options, libs, label = "Codon") {
  ns <- NS(id)
  genomes <- unique(all_exp$organism)
  experiments <- all_exp$name
  tabPanel(
    title = "Codon dwell time", icon = icon("layer-group"),
    sidebarLayout(
      sidebarPanel(
        tabsetPanel(
          tabPanel("Codon",
                   organism_input_select(c("ALL", genomes), ns),
                   experiment_input_select(experiments, ns, browser_options),
                   gene_input_select(ns),
                   tx_input_select(ns),
                   library_input_select(ns, TRUE, libs),
                   fluidRow(column(6, codon_filter_input_select(ns, 1000)),
                            column(6, sliderInput(ns("ratio_thresh"), "Ratio threshold", min = 1.1, max = 3,
                                                  value =1.7, step = 0.2))),
                   codon_score_input_select(ns),
                   checkboxInput(ns("differential"), label = "Differential", value = FALSE),
                   checkboxInput(ns("exclude_start_stop"), label = "Exclude start stop", value = TRUE)
                   )),
        plot_button(ns("go"))
      ),
      mainPanel(
        tabsetPanel(type = "tabs",
                    tabPanel("Codon plot", plotlyOutput(outputId = ns("c"), height = "950px") %>% shinycssloaders::withSpinner(color="#0dc5c1")),
                    tabPanel("Codon table", DTOutput(outputId = ns("codon_table")) %>% shinycssloaders::withSpinner(color="#0dc5c1"))))
  ))
}

codon_server <- function(id, all_experiments, env, df, experiments, tx, cds,
                         libs, org, gene_name_list, rv) {
  moduleServer(
    id,
    function(input, output, session, all_exp = all_experiments) {
      length_table <- reactive(optimizedTranscriptLengths(df(), TRUE, TRUE)) %>%
        bindCache(rv$curval) %>%
        bindEvent(rv$changed, ignoreNULL = TRUE)
      # Update main side panels
      all_is_gene <- TRUE
      study_and_gene_observers(input, output, session)

      # Main plot, this code is only run if 'plot' is pressed
      mainPlotControls <- eventReactive(input$go,
                     click_plot_codon_main_controller(input, tx, cds, libs, df,
                                                      length_table))

      coverage <- reactive(codon_data(mainPlotControls, tx)) %>%
        bindCache(mainPlotControls()$hash_string)
      output$c <- renderPlotly(click_plot_codon_shiny(mainPlotControls(), coverage())) %>%
        bindEvent(coverage(), ignoreInit = FALSE, ignoreNULL = TRUE)
      output$codon_table <- DT::renderDT(coverage(),
                                         extensions = 'Buttons',
                                         filter = "top",
                                         options = list(dom = 'Bfrtip',
                                                        buttons = NULL,
                                                        pageLength = 130))
      return(rv)
    }
  )
}
