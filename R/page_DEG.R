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
                   factor_input_select(ns),
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
        plot_button(ns("go")),
      ),
      mainPanel(
        tabsetPanel(type = "tabs",
                    tabPanel("DEG plot", tagList(uiOutput(outputId = ns("c"), height = "500px") %>% shinycssloaders::withSpinner(color="#0dc5c1"),
                                                 uiOutput(ns("selected_info")))),
                    tabPanel("DEG table", DTOutput(outputId = ns("deg_table")) %>% shinycssloaders::withSpinner(color="#0dc5c1")),
                    tabPanel("GO analysis",
                             fluidRow(
                               column(3, actionButton(ns("run_gorilla"), "Run GORilla", class = "btn btn-primary")),
                               column(9, DTOutput(ns("gorilla_dt")) %>% shinycssloaders::withSpinner(color="#0dc5c1"))))))
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

      factor <- reactive(
        if (nrow(df()) > 1) {
          design <- design(df(), batch.correction.design = TRUE, multi.factor = TRUE)
          names(design) <- stringr::str_to_sentence(design)
          names(design)[design == "rep"] <- "Replicate"
          return(design)
        } else "")
      observeEvent(factor(), factor_update_select(factor))

      cond <- reactive({
        req(isTruthy(input$factor))
        factor <- isolate(input$factor)
        if (nrow(df()) > 1 & factor != "") {
          df()[, input$factor]
        } else ""
      })
      observeEvent(cond(), condition_update_select(cond))

      observeEvent(input$condition, {
        req(isTruthy(input$condition))
        cond_1 <- input$condition[1]
        library_update_select(libs, libs()[cond() == cond_1], "library1")
        }, ignoreNULL = TRUE, ignoreInit = FALSE)

      observeEvent(input$condition, {
        req(isTruthy(input$condition))
        req(length(input$condition) > 1)
        cond_2 <- input$condition[2]
        library_update_select(libs, libs()[cond() == cond_2], "library2")
      }, ignoreNULL = TRUE, ignoreInit = FALSE)

      # Main plot, this code is only run if 'plot' is pressed
      controls <- eventReactive(input$go,
                                click_plot_DEG_main_controller(input, df, libs, factor()))
      model <- reactive(DE_model_from_ctrl(controls)) %>%
        bindCache(controls()$hash_string_pre)
      #
      analysis_dt <- reactive(DE_model_results(model(), controls, gene_name_list())) %>%
        bindCache(controls()$hash_string_full)
      output$c <- renderUI({
        p <- DEG_plot(analysis_dt(), add_search_bar = FALSE,
                      draw_non_regulated = controls()$draw_unregulated,
                      format = controls()$plot_export_format)
        p$x$source <- NS(id)("c")
        event_register(p, "plotly_click")
        DEG_add_search_bar(p)
        }) %>%
        bindCache(controls()$hash_string_plot)


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
                              plot_on_start = TRUE, frames_type = "area", kmer = 1,
                              add_translons = FALSE, host = host)

          google_url <- paste0("https://www.google.com/search?q=", URLencode(selected_key))
          ribocrypt_browser_url <- urls

          tagList(
            strong("Selected ID: "), selected_key, br(), br(),
            strong("Search ID on: "),
            a("RiboCrypt Browser", href = ribocrypt_browser_url, target = "_blank",
              style = "color: blue; text-decoration: underline; margin-right: 10px;"),
            a("Google", href = google_url, target = "_blank",
              style = "color: green; text-decoration: underline; margin-right: 10px;"),
            actionButton(NS(id)("show_description"), "Show Description", icon = icon("info-circle"),
                         class = "btn btn-sm btn-outline-secondary"), br(), br(),
            uiOutput(NS(id)("gene_description")) %>% shinycssloaders::withSpinner(color="#0dc5c1")
          )
        } else if (isolate(input$go) > 0) {
          "Click a point to see details here"
        } else {
          ""
        }
      }) %>% bindEvent(c(input$go, event_data("plotly_click", source = NS(id)("c"))), ignoreInit = TRUE)

      output$deg_table <- DT::renderDT(analysis_dt(),
                                       extensions = 'Buttons',
                                       filter = "top",
                                       options = list(dom = 'Bfrtip',
                                                      buttons = NULL,
                                                      pageLength = 130))

      output$gene_description <- renderUI({
        selected_id <- suppressMessages(event_data("plotly_click", source = NS(id)("c")))
        if (!is.null(selected_id)) {
          selected_key <- selected_id$key
          tx_clicked_raw <- sub("\\(.*", "", selected_key)
          tx_clicked_dt <- gene_name_list()[value == tx_clicked_raw]
          gene_and_symbol <- tx_clicked_dt$label
          gene_subset <- sub("-.*", "", gene_and_symbol)
          format <- ifelse(gene_and_symbol == gene_subset, "ensembl_id", "symbol")
          gene_desc <- ORFik::download_gene_info(gene_subset, organism = organism(df()), by = format)
          if (is.null(gene_desc) || is.na(gene_desc) || gene_desc == "") {
            gene_desc <- "No description available for this gene."
          }
          div(style = "white-space: pre-wrap; font-size: 14px; color: #333; padding-top: 5px;",
              gene_desc)
        }
      }) %>% bindEvent(c(input$show_description, event_data("plotly_click", source = NS(id)("c"))), ignoreInit = TRUE)

      output$gorilla_dt <- renderDT(NULL)

      observeEvent(input$run_gorilla, {
        print("Running gorilla!")
        req(analysis_dt())

        output$gorilla_dt <- renderDT({
          dt <- analysis_dt()
          dt[, external_gene_name := sub("-.*", "", label)]
          gorilla_output_dir <- tempdir()
          gorilla_result <- ORFik:::DEG_gorilla(dt = dt, output_dir = gorilla_output_dir, organism(df()))
          gorilla_result[, url := sprintf('<a href="%s" target="_blank" style="color: blue;">Link</a>', url)]
          datatable(gorilla_result, escape = FALSE, options = list(pageLength = 25))
        })
      })
      return(rv)
    }
  )
}
