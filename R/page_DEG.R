DEG_ui <- function(id, all_exp, browser_options, label = "DEG") {
  ns <- NS(id)
  genomes <- unique(all_exp$organism)
  experiments <- all_exp$name
  tabPanel(
    title = "Differential expression", icon = icon("layer-group"),
    sidebarLayout(
      sidebarPanel(
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
      ),
      mainPanel(
        jqui_resizable(plotlyOutput(outputId = ns("c"), height = "500px")) %>% shinycssloaders::withSpinner(color="#0dc5c1"),
        uiOutput(ns("selected_info"))
        )
    )
  )
}


DEG_server <- function(id, all_experiments, env, df, experiments, libs,
                       org, gene_name_list, rv) {
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
        }, ignoreNULL = TRUE, ignoreInit = FALSE)

      observeEvent(cond(), {
        cond_2 <- unique(cond())
        req(length(cond_2) > 1)
        cond_2 <- cond_2[2]
        library_update_select(libs, libs()[cond() == cond_2], "library2")
      }, ignoreNULL = TRUE, ignoreInit = FALSE)

      # Main plot, this code is only run if 'plot' is pressed
      controls <- eventReactive(input$go,
                                click_plot_DEG_main_controller(input, df, libs))
      model <- reactive(DE_model_from_ctrl(controls)) %>%
        bindCache(controls()$hash_string_pre)
      #
      analysis_dt <- reactive(DE_model_results(model(), controls, gene_name_list())) %>%
        bindCache(controls()$hash_string_full)
      output$c <- renderPlotly({
        p <- DEG_plot(analysis_dt(), draw_non_regulated = controls()$draw_unregulated)
        p$x$source <- NS(id)("c")
        event_register(p, "plotly_click")
        p
        }) %>%
        bindCache(controls()$hash_string_full)


      output$selected_info <- renderUI({
        message("Rendering DEG text UI")
        selected_id <- suppressMessages(event_data("plotly_click", source = NS(id)("c")))
        if (!is.null(selected_id)) {
          selected_key <- selected_id$key
          tx_clicked_raw <- sub("\\(.*", "", selected_key)
          tx_clicked_dt <- gene_name_list()[value == tx_clicked_raw]
          tx_clicked <- tx_clicked_dt$value


          host <- getHostFromURL(session)
          exp <- name(df())
          gene_and_symbol <- tx_clicked_dt$label
          libraries_to_use <- c(input$library1, input$library2)
          urls <- make_rc_url(symbol = NULL, gene_id = gene_and_symbol, tx_id = tx_clicked,
                              exp = exp,
                              libraries = libraries_to_use, leader_extension = 0, trailer_extension = 0,
                              viewMode = FALSE, other_tx = FALSE,
                              plot_on_start = TRUE, frames_type = "area", kmer=1,
                              add_translons = FALSE, host = host)
          google_url <- paste0("https://www.google.com/search?q=", URLencode(selected_key))
          ribocrypt_browser_url <- urls

          # Create two clickable hyperlinks
          tagList(
            strong("Selected ID: "), selected_key, br(), br(),
            strong("Search ID on: "),
            a("RiboCrypt Browser", href = ribocrypt_browser_url, target = "_blank", style = "color: blue; text-decoration: underline; margin-right: 10px;"),
            a("Google", href = google_url, target = "_blank", style = "color: green; text-decoration: underline;")
          )
        } else if (isolate(input$go) > 0) {
          "Click a point to see details here"
        } else {""}
      }) %>% bindEvent(c(input$go, event_data("plotly_click", source = NS(id)("c"))), ignoreInit = TRUE)


      return(rv)
    }
  )
}
