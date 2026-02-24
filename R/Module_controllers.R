### NGLVieweR (protein structures) ###
# TODO: Move as much as possible of protein stuff out of page_browser

module_protein <- function(input, output, gene_name_list, session) {
  with(rlang::caller_env(), {
    # Setup reactive values needed for structure viewer
    dynamicVisible <- reactiveVal(FALSE)
    selectedRegion <- reactiveVal(NULL)
    selectedTX <- reactiveVal(NULL)
    if (FALSE) cds <- NULL # Avoid biocCheck error

    # When user clicks on region
    # start displaying structure viewer
    # and set selected structure to one which was clicked
    observeEvent(input$selectedRegion, {
      req(input$selectedRegion, input$useCustomRegions)
      selectedRegion(input$selectedRegion)
      selectedTX(input$tx)
      dynamicVisible(TRUE)
    }, priority = 1)

    # When user clicks close button
    # stop displaying structure viewer
    # and set selected structure to NULL
    observeEvent(input$dynamicClose, {
      selectedRegion(NULL)
      dynamicVisible(FALSE)
    })

    # NGL viewer widget
    protein_structure_dir <- reactive({
      file.path(refFolder(df()), "protein_structure_predictions")
    })
    region_dir <- reactive({
      file.path(protein_structure_dir(), selectedTX())
    })


    structure_variants <- reactive({
      req(selectedRegion(), selectedRegion() != "...")
      get_protein_structure(gene_name_list(), selectedRegion(), selectedTX(), region_dir(),
                            protein_structure_dir(), df(), mainPlotControls()$customRegions)
    }) %>% bindCache(selectedRegion(), selectedTX) %>%
      bindEvent(selectedRegion(), ignoreNULL = TRUE, ignoreInit = FALSE)

    # Get the Ribo-seq profile (we select first library for now)
    selectedRegionProfile <- reactive({
      req(selectedRegion(), selectedRegion() != "...", input$useCustomRegions,
          structure_variants())
      getCoverageProtein(mainPlotControls()$reads[[1]], cds(), mainPlotControls()$customRegions,
                         selectedRegion(), input$log_scale_protein, id)

    }) %>% bindCache(selectedRegion(), selectedTX())%>%
      bindEvent(structure_variants(), ignoreNULL = TRUE, ignoreInit = FALSE)

    # Variable UI logic
    output$dynamic <- renderNGLVieweR(protein_struct_render(selectedRegionProfile(),
                                                            attr(structure_variants(), "selected")))
    # output$dynamic <- renderR3dmol(protein_struct_render(selectedRegionProfile, selected_variant))
    output$proteinStruct <- renderUI(
      protein_struct_plot(
        selectedRegion,
        selectedRegionProfile,
        dynamicVisible,
        session,
        structure_variants
      )
    )
  })
}

module_additional_browser <- function(input, output, session) {
  with(rlang::caller_env(), {
    # Protein display
    module_protein(input, output, gene_name_list, session)

    output$download_plot_html <- downloadHandler(
      filename = function() {
        paste0("RiboCrypt_", isolate(input$tx), Sys.Date(), ".html")
      },
      content = function(file) {
        htmlwidgets::saveWidget(as_widget(browser_plot()), file)
      }
    )

    observeEvent(input$toggle_settings, {
      # Toggle visibility by adding/removing 'hidden' class
      shinyjs::toggleClass(id = "floating_settings", class = "hidden")
    })

    observeEvent(input$myInput_copy, {
      showNotification(paste("Copied", nchar(input$myInput_copy), "nt to clipboard"), type = "message")
    })

    observe({
      if(!isTruthy(input$go)) {
        shinyjs::hideElement(id = "download_plot_html")
      } else {
        shinyjs::showElement(id = "download_plot_html")
      }
    })

    output$d <- renderPlotly({
      req(input$expression_plot == TRUE)
      click_plot_boxplot(mainPlotControls, session)}) %>%
      bindCache(mainPlotControls()$hash_expression)
  })
}

module_additional_megabrowser <- function(input, output, session) {
  with(rlang::caller_env(), {
    selected_enrich_filters <- reactiveVal(NULL)
    observeEvent(meta_and_clusters(), {
      selected_enrich_filters(NULL)
    }, ignoreInit = TRUE)

    observeEvent(plotly::event_data("plotly_click", source = "mb_enrich"), {
      ed <- plotly::event_data("plotly_click", source = "mb_enrich")
      req(!is.null(ed))

      category <- as.character(ed$x)
      cluster <- NULL
      if (!is.null(ed$customdata) && !identical(ed$customdata, "counts")) {
        cluster <- as.character(ed$customdata)
      }

      selected_enrich_filters(list(category = category, cluster = cluster))
      updateTabsetPanel(session, "mb_tabs", selected = "Result table")
      shinyjs::runjs(
        sprintf(
          "document.getElementById('%s').scrollIntoView({behavior:'smooth'});",
          ns("result_table")
        )
      )
    })

    observeEvent(plotly::event_data("plotly_clickannotation", source = "mb_sidebar"), {
      ed <- plotly::event_data("plotly_clickannotation", source = "mb_sidebar")
      req(!is.null(ed))

      cluster_val <- NULL
      if (is.list(ed) && !is.null(ed$text)) {
        cluster_val <- as.character(ed$text)
      } else if (is.data.frame(ed) && "text" %in% names(ed)) {
        cluster_val <- as.character(ed$text[[1]])
      } else if (is.list(ed) && !is.null(ed$annotation) && !is.null(ed$annotation$text)) {
        cluster_val <- as.character(ed$annotation$text)
      } else if (is.data.frame(ed) && "annotation.text" %in% names(ed)) {
        cluster_val <- as.character(ed[["annotation.text"]][[1]])
      }
      req(!is.null(cluster_val))

      selected_enrich_filters(list(category = NULL, cluster = cluster_val))
      updateTabsetPanel(session, "mb_tabs", selected = "Result table")
      shinyjs::runjs(
        sprintf(
          "document.getElementById('%s').scrollIntoView({behavior:'smooth'});",
          ns("result_table")
        )
      )
    })

    selected_cluster_from_filtered <- reactive({
      filt <- selected_enrich_filters()
      if (is.null(filt)) return(NULL)
      if (!is.null(filt$cluster)) return(as.character(filt$cluster))

      tbl <- filtered_meta_table()
      if ("cluster" %in% names(tbl)) {
        vals <- unique(tbl$cluster)
      } else if ("Cluster" %in% names(tbl)) {
        vals <- unique(tbl$Cluster)
      } else {
        vals <- character(0)
      }
      if (length(vals) == 1) return(as.character(vals))
      NULL
    })

    output$result_table_controls <- renderUI({
      filt <- selected_enrich_filters()
      if (is.null(filt)) return(NULL)

      show_cluster_btn <- FALSE
      cluster_val <- selected_cluster_from_filtered()
      if (!is.null(cluster_val)) {
        already_full_cluster <- is.null(filt$category) && !is.null(filt$cluster)
        show_cluster_btn <- !already_full_cluster
      }

      tagList(
        actionButton(ns("reset_result_table"), "Show full table"),
        if (show_cluster_btn) actionButton(ns("expand_cluster"), "Show full cluster") else NULL
      )
    })

    observeEvent(input$reset_result_table, {
      selected_enrich_filters(NULL)
    }, ignoreInit = TRUE)

    observeEvent(input$expand_cluster, {
      cluster_val <- selected_cluster_from_filtered()
      req(!is.null(cluster_val))
      selected_enrich_filters(list(category = NULL, cluster = cluster_val))
    }, ignoreInit = TRUE)

    filtered_meta_table <- reactive({
      tbl <- allsamples_meta_table(meta_and_clusters())
      filt <- selected_enrich_filters()
      if (is.null(filt)) return(tbl)
      if (!is.null(filt$category) && "grouping" %in% names(tbl)) {
        tbl <- tbl[as.character(grouping) %in% filt$category]
      }
      if (!is.null(filt$cluster)) {
        if ("cluster" %in% names(tbl)) {
          tbl <- tbl[as.character(cluster) %in% filt$cluster]
        } else if ("Cluster" %in% names(tbl)) {
          tbl <- tbl[as.character(Cluster) %in% filt$cluster]
        }
      }
      tbl
    })

    output$result_table <- renderDT(filtered_meta_table(),
                                    extensions = 'Buttons', filter = "top",
                                    options = list(dom = 'Bfrtip',
                                                   buttons = NULL,
                                                   pageLength = 130)) %>%
      bindEvent(meta_and_clusters(), selected_enrich_filters(),
                ignoreInit = FALSE,
                ignoreNULL = TRUE)
  })
}

#' This function sets up default backend for genome specific reactives
#'
#' It is a rlang module for all submodules.\cr
#' It has branch points for setup to be more flexible for modules.
#' @noRd
study_and_gene_observers <- function(input, output, session) {
  with(rlang::caller_env(), {
    # Checks for which flags to set from parent function
    if (!exists("all_is_gene", mode = "logical")) all_is_gene <- FALSE
    if (!exists("uses_gene", mode = "logical")) uses_gene <- TRUE
    if (!exists("uses_libs", mode = "logical")) uses_libs <- TRUE
    if (!exists("env")) env <- new.env()

    observe(if (rv$genome != input$genome & isTruthy(input$genome)) {
      message("Changing org from org switch in other page")
      rv$genome <- input$genome},
      priority = 2) %>%
      bindEvent(input$genome, ignoreInit = TRUE, ignoreNULL = TRUE)
    observe(if (rv$exp != input$dff && input$dff != "") {
      message("Setting rv from page: ", id)
      rv$exp <- input$dff}) %>%
      bindEvent(input$dff, ignoreInit = TRUE, ignoreNULL = TRUE)

    observeEvent(rv$exp, if (rv$exp != input$dff) {
      experiment_update_select(org, all_exp, experiments, rv$exp)},
      ignoreInit = TRUE, ignoreNULL = TRUE)

    observe(if (rv$genome != input$genome) {
      updateSelectizeInput(
        inputId = "genome",
        choices = c("ALL", unique(all_exp$organism)),
        selected = rv$genome,
        server = TRUE
      )}, priority = 1) %>%
      bindEvent(rv$genome, ignoreInit = TRUE, ignoreNULL = TRUE)
    experiment_update_select_isolated(isolate(org()), all_exp, experiments,
                                      isolate(rv$exp))

    observeEvent(org(), if (org() != input$genome & isTruthy(input$genome)) {
      message("Changing exp from org switch")
      experiment_update_select(org, all_exp, experiments)},
      ignoreInit = TRUE, ignoreNULL = TRUE)
    # Gene & tx updaters
    if (all_is_gene) {
      updateSelectizeInput(
        inputId = "gene",
        choices = c("all", unique(isolate(gene_name_list())[,2][[1]])),
        selected = "all",
        server = TRUE
      )
      observeEvent(gene_name_list(), gene_update_select_heatmap(gene_name_list),
                   ignoreInit = FALSE)
      observeEvent(input$gene, {
        req(input$gene != "")
        tx_update_select(isolate(input$gene), gene_name_list, "all", page = id)
      }, ignoreNULL = TRUE, ignoreInit = FALSE)

    } else if (uses_gene) {
      print(id)
      choices <- unique(isolate(gene_name_list())[,2][[1]])
      # Init round gene
      if (id == "browser_allsamp") {
        print("Updating metabrowser gene set")
        gene_update_select_internal(isolate(gene_name_list()), selected = browser_options["default_gene_meta"])
        gene_update_select_internal(NULL, choices = c("", choices),
                                    id = "other_gene")
        observeEvent(gene_name_list(), gene_update_select(gene_name_list, "",
                                                          id = "other_gene"),
                     ignoreNULL = TRUE, ignoreInit = FALSE, priority = 6)
      } else {
        gene_update_select_internal(isolate(gene_name_list()), selected = browser_options["default_gene"])
      }
      # Non init round gene
      observeEvent(gene_name_list(), {
        selected <- unique(isolate(gene_name_list())[,2][[1]])[1]
        gene_update_select(gene_name_list, selected = selected)
      }, ignoreNULL = TRUE, ignoreInit = TRUE, priority = 5)

      # Tx id update
      # Init round
      if (id == "browser_allsamp") {
        tx_update_select_isolated(browser_options["default_gene_meta"], isolate(gene_name_list()),
                                  selected = browser_options["default_isoform_meta"], page = id)
      } else {
        tx_update_select_isolated(browser_options["default_gene"], isolate(gene_name_list()),
                                  selected = browser_options["default_isoform"], page = id)
      }
      # Non int rounds
      observeEvent(input$gene, {
        req(input$gene != "")
        if (id != "browser_allsamp") {
          req(!(input$tx %in% c("",
                isolate(gene_name_list())[label == input$gene,]$value)))
        }
        print(paste("Page:", id, "(General observer)"))
        tx_update_select_isolated(isolate(input$gene), isolate(gene_name_list()), page = id)
        },
        ignoreNULL = TRUE, ignoreInit = TRUE, priority = -15)
    }
    if (uses_libs) {
      if (!exists("init_round") && exists("browser_options")) {
        selected_libs <- libraries_string_split(browser_options["default_libs"], isolate(libs()))
        library_update_select_safe(isolate(libs()), selected_libs)
      }
      observeEvent(libs(), library_update_select(libs),
                   ignoreNULL = TRUE, ignoreInit = TRUE)
      if (id == "browser") {
        observeEvent(input$select_all_btn, {
          print("Pressed select all libs")
          library_update_select(libs, selected = libs())
        }, ignoreNULL = TRUE, ignoreInit = TRUE)  # Ensure the event is triggered on every click, including after initialization
      } else if (id == "Codon") {
        observeEvent(libs(), library_update_select(libs, id = "background"),
                     ignoreNULL = TRUE, ignoreInit = TRUE)
      }
    }
    browser_specific_url_checker()

    init_round <- FALSE
  }
  )
}

org_and_study_changed_checker <- function(input, output, session) {
  with(rlang::caller_env(), {
    ## Static values
    experiments <- all_exp$name
    ## Set reactive values
    org <- reactiveVal("ALL")

    df <- reactiveVal(get_exp(browser_options["default_experiment"],
                              experiments, without_readlengths_env, exps_dir))
    df_with <- reactiveVal(get_exp(browser_options["default_experiment"],
                              experiments, with_readlengths_env, exps_dir))
    if (nrow(all_exp_meta) > 0) {
      df_meta <- reactiveVal({
        if(name(exp_init_meta) == browser_options["default_experiment_meta"]) {
          print(paste("Loading exp:", name(exp_init_meta)))
          print("- Init experiment loaded")
          exp_init_meta
        } else {get_exp(browser_options["default_experiment_meta"],
                        all_exp_meta$name, .GlobalEnv, exps_dir)}
      })
    } else print("No MegaBrowser exps given, ignoring MegaBrowser exp.")

    libs <- reactive(bamVarName(df()))
    # The shared reactive values (rv)
    # This must be passed to all submodules that share experiment
    rv <- reactiveValues(lstval=isolate(df())@txdb,
                         curval=isolate(df())@txdb,
                         initval=isolate(df())@txdb,
                         genome = "ALL",
                         exp = browser_options["default_experiment"],
                         changed=isolate(df())@txdb != exp_init@txdb)
    observe(update_rv_changed(rv), priority = 1) %>%
      bindEvent(rv$curval, ignoreInit = TRUE)
    observe({update_rv(rv, df)}) %>%
      bindEvent(df(), ignoreInit = TRUE)
    observe({update_rv(rv, df_with)}) %>%
      bindEvent(df_with(), ignoreInit = TRUE)

    observe(if (org() != rv$genome) org(rv$genome)) %>%
      bindEvent(rv$genome, ignoreInit = TRUE, ignoreNULL = TRUE)
    observe({df(get_exp(rv$exp, experiments, without_readlengths_env, exps_dir))}) %>%
      bindEvent(rv$exp, ignoreInit = TRUE, ignoreNULL = TRUE)
    observe({df_with(get_exp(rv$exp, experiments, with_readlengths_env, exps_dir))}) %>%
      bindEvent(rv$exp, ignoreInit = TRUE, ignoreNULL = TRUE)

    # Annotation change reactives
    tx <- reactive({
      if(rv$curval == rv$initval) {message("Settings tx to init:"); tx_init}
      else {loadRegion(isolate(df()))}}) %>%
      bindCache(rv$curval) %>%
      bindEvent(rv$changed, ignoreNULL = TRUE)
    cds <- reactive({
      if(rv$curval == rv$initval) {
        message("Settings cds to init:"); cds_init}
      else {loadRegion(isolate(df()), "cds")}}) %>%
      bindCache(rv$curval) %>%
      bindEvent(rv$changed, ignoreNULL = TRUE)
    # gene_name_list <- reactiveVal(names_init)
    gene_name_list <- reactive({
      if(rv$curval == rv$initval) {
        message("Settings gene_list to init:"); names_init}
      else {get_gene_name_categories(df())}}) %>%
      bindCache(rv$curval) %>%
      bindEvent(rv$changed, ignoreNULL = TRUE)

    cat("Pre server modules: "); print(round(Sys.time() - time_before, 2))
  }
  )
}

allsamples_observer_controller <- function(input, output, session) {
  with(rlang::caller_env(), {
  org <- reactiveVal("ALL")

  rv <- reactiveValues(lstval=isolate(df())@txdb,
                       curval=isolate(df())@txdb,
                       genome = "ALL",
                       exp = name(isolate(df())),
                       changed=isolate(df())@txdb != exp_init@txdb)
  observe(update_rv_changed(rv), priority = 1) %>%
    bindEvent(rv$curval, ignoreInit = TRUE)
  observe({update_rv(rv, df)}) %>%
    bindEvent(df(), ignoreInit = TRUE)

  observe(if (org() != rv$genome) org(rv$genome)) %>%
    bindEvent(rv$genome, ignoreInit = TRUE, ignoreNULL = TRUE)
  observe({df(get_exp(rv$exp, experiments, envExp(df()), exps_dir))}) %>%
    bindEvent(rv$exp, ignoreInit = TRUE, ignoreNULL = TRUE)

  gene_name_list <- reactive({
    if(rv$changed == FALSE) {names_init}
    else {get_gene_name_categories(df())}}) %>%
    bindCache(rv$curval) %>%
    bindEvent(rv$changed)
  motif_name_list <- reactive({
    names(meta_motif_files(df()))}) %>%
    bindCache(rv$curval) %>%
    bindEvent(rv$changed)

  init_motfis <- names(meta_motif_files(isolate(df())))
  motif_update_select(init_motfis)
  observeEvent(motif_name_list(), motif_update_select(motif_name_list()),
               ignoreInit = TRUE)
  uses_libs <- FALSE # Assign for line below
  study_and_gene_observers(input, output, session)
  })
}

meta_motif_files <- function(df) {
  motif_files <- list.files(file.path(resFolder(df), "meta_collection_tables"), pattern = "fst$", full.names = TRUE)
  names(motif_files) <- sub("\\.fst$", "", gsub(".*_", "", basename(motif_files)))
  return(motif_files)
}
