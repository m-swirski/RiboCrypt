### NGLVieweR (protein structures) ###
# TODO: Move as much as possible of protein stuff out of page_browser
module_protein <- function(input, output, gene_name_list, session) {
  with(rlang::caller_env(), {
    # Setup reactive values needed for structure viewer
    dynamicVisible <- reactiveVal(FALSE)
    selectedRegion <- reactiveVal(NULL)
    if (FALSE) cds <- NULL # Avoid biocCheck error
    # Get the Ribo-seq prfile (we select first library for now)
    selectedRegionProfile <- reactive({
      req(selectedRegion())
      result <- cds()[names(cds()) == selectedRegion()] %>%
        getRiboProfile(mainPlotControls()$reads[[1]]) %>%
        (function (x) {
          x$count[seq.int(1, length(x$count), 3)]
        })()
    })

    # When user clicks on region
    # start displaying structure viewer
    # and set selected structure to one which was clicked
    observeEvent(input$selectedRegion, {
      req(input$selectedRegion)
      selectedRegion(input$selectedRegion)
      dynamicVisible(TRUE)
    })
    # When user clicks close button
    # stop displaying structure viewer
    # and set selected structure to NULL
    observeEvent(input$dynamicClose, {
      selectedRegion(NULL)
      dynamicVisible(FALSE)
    })
    # Setup 3dbeacons model download
    observeEvent(beacons_structures(), {
      print("Fetching structures..")
      tmp_paths <- unname(beacons_structures())
      names(tmp_paths) <- beacons_results()
      mapply(
        function(x) {
          httr::GET(x, httr::write_disk(tmp_paths[x], overwrite = TRUE))
        },
        beacons_results()
      )
    })

    # NGL viewer widget
    protein_structure_dir <- reactive({
      file.path(dirname(df()@fafile), "protein_structure_predictions")
    })
    region_dir <- reactive({
      file.path(protein_structure_dir(), selectedRegion())
    })
    on_disk_structures <- reactive({
      paths <- file.path(region_dir(), list.files(region_dir()))

      path_labels <- mapply(
        function(x) {
          str_sub(x, start = gregexpr("/", x) %>% unlist() %>% last() + 1)
        },
        list.files(region_dir())
      )

      result <- paths
      names(result) <- path_labels
      result
    })
    uniprot_id <- reactive(
      gene_name_list()[gene_name_list()$value == selectedRegion()]$uniprot_id
    )

    beacons_results <- reactive({
      print("Fetching structures URLs..")
      model_urls <-
        fetch_summary(uniprot_id()) %>%
        model_urls_from_summary()
      # Convert alphafold cif to pdb (NGLviewR does not work with cif there)
      alphafold_id <- grep("https://alphafold.ebi.ac.uk/files/", model_urls)
      if (length(alphafold_id) > 0) {
        model_urls[alphafold_id] <- gsub("cif$", "pdb", model_urls[alphafold_id])
      }
      model_urls
    }) %>%
      bindCache(uniprot_id()) %>%
      bindEvent(uniprot_id(), ignoreNULL = TRUE, ignoreInit = TRUE)

    beacons_structures <- reactive({
      req(beacons_results())
      results <- beacons_results()

      model_labels <- mapply(
        function(x) {
          str_sub(x, start = gregexpr("/", x) %>% unlist() %>% last() + 1)
        },
        results
      )
      model_paths <- rep(tempfile(pattern = "structure"), length(model_labels))
      model_paths <- paste0(model_paths, ".", file_ext(results))
      result <- model_paths
      names(result) <- model_labels
      result
    })
    structure_variants <- reactive({
      print("Structures fetched")
      print(beacons_structures())
      append(on_disk_structures(), beacons_structures())
    })
    selected_variant <- reactive({
      req(structure_variants())
      if (is.null(input$structureViewerSelector)) {
        head(structure_variants())
      } else input$structureViewerSelector
    })
    # Variable UI logic
    output$dynamic <- renderNGLVieweR(protein_struct_render(selectedRegionProfile, selected_variant))
    # output$dynamic <- renderR3dmol(protein_struct_render(selectedRegionProfile, selected_variant))
    output$variableUi <- renderUI(
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

### NGLVieweR (protein structures) ###
# TODO: Move as much as possible of protein stuff out of page_browser
### NGLVieweR (protein structures) ###
# TODO: Move as much as possible of protein stuff out of page_browser
study_and_gene_observers <- function(input, output, session) {
  with(rlang::caller_env(), {
    if (!exists("all_is_gene", mode = "logical")) all_is_gene <- FALSE
    if (!exists("uses_gene", mode = "logical")) uses_gene <- TRUE
    if (!exists("uses_libs", mode = "logical")) uses_libs <- TRUE
    if (!exists("env")) env <- new.env()

    observe(if (rv$genome != input$genome & input$genome != "") {
      rv$genome <- input$genome},
      priority = 2) %>%
      bindEvent(input$genome, ignoreInit = TRUE, ignoreNULL = TRUE)
    observe(if (rv$exp != input$dff & input$dff != "") rv$exp <- input$dff) %>%
      bindEvent(input$dff, ignoreInit = TRUE, ignoreNULL = TRUE)
    observe(if (rv$genome != input$genome) {
      updateSelectizeInput(
        inputId = "genome",
        choices = c("ALL", unique(all_exp$organism)),
        selected = rv$genome,
        server = TRUE
      )}, priority = 1) %>%
      bindEvent(rv$genome, ignoreInit = TRUE, ignoreNULL = TRUE)

    observeEvent(rv$exp, if (rv$exp != input$dff) {
      experiment_update_select(org, all_exp, experiments, rv$exp)},
      ignoreInit = TRUE, ignoreNULL = TRUE)

    observeEvent(org(), if (org() != input$genome & input$genome != "") {
      experiment_update_select(org, all_exp, experiments)},
      ignoreInit = TRUE, ignoreNULL = TRUE)
    if (all_is_gene) {
      updateSelectizeInput(
        inputId = "gene",
        choices = c("all", unique(isolate(gene_name_list())[,2][[1]])),
        selected = "all",
        server = TRUE
      )
      observeEvent(gene_name_list(), gene_update_select_heatmap(gene_name_list),
                   ignoreInit = TRUE)
      observeEvent(input$gene, {
        req(input$gene != "")
        tx_update_select(isolate(input$gene), gene_name_list, "all")
      }, ignoreNULL = TRUE, ignoreInit = TRUE)

    } else if (uses_gene) {
      print(id)
      if (id == "browser_allsamp") {
        print("Updating metabrowser gene set")
        choices <- unique(isolate(gene_name_list())[,2][[1]])
        updateSelectizeInput(
          inputId = "gene",
          choices = choices,
          selected = choices[1],
          server = TRUE
        )
      }
      # TODO: decide if updateSelectizeInput should be on top here or not
      observeEvent(gene_name_list(), gene_update_select(gene_name_list),
                   ignoreNULL = TRUE, ignoreInit = TRUE, priority = 5)

      check_url_for_basic_parameters()
      observeEvent(input$gene, {
        req(input$gene != "")
        if (id != "browser_allsamp") {
          req(!(input$tx %in% c("",
                isolate(gene_name_list())[label == input$gene,]$value)))
        }

        print(paste("Page:", id, "(General observer)"))
        tx_update_select(isolate(input$gene), gene_name_list)},
        ignoreNULL = TRUE, ignoreInit = TRUE, priority = -15)
      print("Updating browser gene set")
      choices <- unique(isolate(gene_name_list())[,2][[1]])
      updateSelectizeInput(
        inputId = "gene",
        choices = choices,
        selected = isolate(input$gene),
        server = TRUE
      )

    }

    if (uses_libs) {
      observeEvent(libs(), library_update_select(libs),
                   ignoreNULL = TRUE, ignoreInit = FALSE)
    }
    check_url_for_go_on_init()
    init_round <- FALSE
  }
  )
}

org_and_study_changed_checker <- function(input, output, session) {
  with(rlang::caller_env(), {
    cat("Server startup: "); print(round(Sys.time() - time_before, 2))
    # browser()
    ## Static values
    experiments <- all_exp$name
    ## Set reactive values
    org <- reactiveVal("ALL")
    # Store current and last genome
    df <- reactiveVal(get_exp(browser_options["default_experiment"],
                              experiments, without_readlengths_env))
    df_with <- reactiveVal(get_exp(browser_options["default_experiment"],
                              experiments, with_readlengths_env))
    df_meta <- reactiveVal(get_exp(browser_options["default_experiment_meta"],
                                   all_exp_meta$name, .GlobalEnv))
    libs <- reactive(bamVarName(df()))
    # The shared reactive values (rv)
    # This must be passed to all submodules
    rv <- reactiveValues(lstval=isolate(df())@txdb,
                         curval=isolate(df())@txdb,
                         genome = "ALL",
                         exp = browser_options["default_experiment"],
                         changed=FALSE)
    # Annotation change reactives
    tx <- reactive(loadRegion(isolate(df()))) %>%
      bindCache(rv$curval) %>%
      bindEvent(rv$changed, ignoreNULL = TRUE)
    cds <- reactive(loadRegion(isolate(df()), "cds")) %>%
      bindCache(rv$curval) %>%
      bindEvent(rv$changed, ignoreNULL = TRUE)
    gene_name_list <- reactiveVal(names_init)
    init_round <- TRUE
    gene_name_list <- reactive({if (init_round) {names_init}
      else get_gene_name_categories(df())}) %>%
      bindCache(rv$curval) %>%
      bindEvent(rv$changed, ignoreNULL = TRUE)
    init_round <- FALSE
    # observe(gene_name_list()) %>%
    #   bindEvent(rv$changed, ignoreInit = TRUE)
    # Observers
    observe(update_rv_changed(rv), priority = 1) %>%
      bindEvent(rv$curval, ignoreInit = TRUE)
    observe({update_rv(rv, df)}) %>%
      bindEvent(df(), ignoreInit = TRUE)
    observe({update_rv(rv, df_with)}) %>%
      bindEvent(df_with(), ignoreInit = TRUE)

    observe(if (org() != rv$genome) org(rv$genome)) %>%
      bindEvent(rv$genome, ignoreInit = TRUE, ignoreNULL = TRUE)
    observe({df(get_exp(rv$exp, experiments, without_readlengths_env))}) %>%
      bindEvent(rv$exp, ignoreInit = TRUE, ignoreNULL = TRUE)
    observe({df_with(get_exp(rv$exp, experiments, with_readlengths_env))}) %>%
      bindEvent(rv$exp, ignoreInit = TRUE, ignoreNULL = TRUE)

    cat("Pre modules: "); print(round(Sys.time() - time_before, 2))
  }
  )
}



