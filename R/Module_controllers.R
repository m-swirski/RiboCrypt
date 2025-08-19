### NGLVieweR (protein structures) ###
# TODO: Move as much as possible of protein stuff out of page_browser

module_protein <- function(input, output, gene_name_list, session) {
  with(rlang::caller_env(), {
    # Setup reactive values needed for structure viewer
    dynamicVisible <- reactiveVal(FALSE)
    selectedRegion <- reactiveVal(NULL)
    selectedTX <- reactiveVal(NULL)
    if (FALSE) cds <- NULL # Avoid biocCheck error

    # Get the Ribo-seq prfile (we select first library for now)
    selectedRegionProfile <- reactive({
      req(selectedRegion(), input$useCustomRegions)
      req(selectedRegion() != "...")
      if (id != "browser") {
        return(rep(1, 1e6))
      }
      coverage_region <- NULL
      uorf_clicked <- length(grep("U[0-9]+$", input$selectedRegion)) == 1
      translon_clicked <- length(grep("T[0-9]+$", input$selectedRegion)) == 1
      if (uorf_clicked) {
        print("- Searching for local uorf protein structure:")
        print(input$selectedRegion)
        print(paste("In tx:", input$tx))
      }
      if (translon_clicked) {
        print("- Searching for local translon protein structure:")
        print(input$selectedRegion)
        print(paste("In tx:", input$tx))
      }
      coverage_region <- c(cds(), mainPlotControls()$customRegions)
      coverage_region <- coverage_region[names(coverage_region) == selectedRegion()]
      stopifnot(length(coverage_region) == 1)

      result <- coverage_region %>%
        getRiboProfile(mainPlotControls()$reads[[1]]) %>%
        (function (x) {
          x$count[seq.int(1, length(x$count), 3)]
        })()
      if (input$log_scale_protein) {
        result <- floor(log2(result))
        result[!is.finite(result)] <- 0
      }
      result
    })

    # When user clicks on region
    # start displaying structure viewer
    # and set selected structure to one which was clicked
    observeEvent(input$selectedRegion, {
      req(input$selectedRegion, input$useCustomRegions)
      selectedRegion(input$selectedRegion)
      selectedTX(input$tx)
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
      file.path(refFolder(df()), "protein_structure_predictions")
    })
    region_dir <- reactive({
      file.path(protein_structure_dir(), selectedTX())
    })
    on_disk_structures <- reactive({
      paths <- character()
      if (!isTruthy(selectedRegion())) return(paths)
      pdb_input <- grepl("\\.pdb$", selectedRegion())
      uorf_clicked <- length(grep("^U[0-9]+$", input$selectedRegion)) == 1
      translon_clicked <- length(grep("^T[0-9]+$", input$selectedRegion)) == 1
      if (pdb_input) {
        paths <- selectedRegion()
      } else if (uorf_clicked) {
        paths <- file.path(region_dir(), list.files(region_dir()))
        paths <- paths[grep(paste0("^uorf_",selectedRegion(), ".pdb$"), basename(paths))]
        if (length(paths) == 0) warning("No local protein structure for this uORF!")
      } else if (translon_clicked) {
        linker_file <- file.path(refFolder(df()), "predicted_translons",
                                 "predicted_translons_with_sequence_pep_linker.fst")
        if (!file.exists(linker_file)) {
          print("No translon peptide linker file for organism!")
          return("")
        }
        selected_grl <- mainPlotControls()$customRegions[]
        selected_as_coord <- as.character(selected_grl[names(selected_grl) == selectedRegion()])
        linker_dt <- fst::read_fst(linker_file, as.data.table = TRUE)

        paths <- coordinates_to_pep_id_path(selected_as_coord, linker_dt, protein_structure_dir())
      } else {
        paths <- paths[-grep("uorf", paths)] # Remove uorf structures
      }


      path_labels <- mapply(
        function(x) {
          str_sub(x, start = gregexpr("/", x) %>% unlist() %>% last() + 1)
        },
        basename(paths)
      )
      result <- paths
      names(result) <- path_labels
      result
    })
    uniprot_id <- reactive({
      gene_name_list()[gene_name_list()$value == selectedRegion()]$uniprot_id
    }

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
      print("Selecting local or online pdb")
      online <- !is(try(beacons_results(), silent = TRUE), "try-error")
      if (online && isTruthy(beacons_results())) {
        print(beacons_structures())
        res <- beacons_structures()
        print("Structures fetched online")
      } else {
        res <- on_disk_structures()
        print("Structures fetched local")
      }
      return(res)
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

    observeEvent(TRUE, {
      experiment_update_select(org, all_exp, experiments, rv$exp)
    }, once = TRUE)

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
        tx_update_select(isolate(input$gene), gene_name_list, "all", page = id)
      }, ignoreNULL = TRUE, ignoreInit = TRUE)

    } else if (uses_gene) {
      print(id)
      choices <- unique(isolate(gene_name_list())[,2][[1]])
      if (id == "browser_allsamp") {
        print("Updating metabrowser gene set")

        gene_update_select_internal(NULL, choices = choices,
                                    id = "gene")
        gene_update_select_internal(NULL, choices = c("", choices),
                                    id = "other_gene")
      }
      # TODO: decide if updateSelectizeInput should be on top here or not
      observeEvent(gene_name_list(), gene_update_select(gene_name_list),
                   ignoreNULL = TRUE, ignoreInit = TRUE, priority = 5)
      observeEvent(gene_name_list(), gene_update_select(gene_name_list, "",
                                                        id = "other_gene"),
                   ignoreNULL = TRUE, ignoreInit = TRUE, priority = 6)

      check_url_for_basic_parameters()
      observeEvent(input$gene, {
        req(input$gene != "")
        if (id != "browser_allsamp") {
          req(!(input$tx %in% c("",
                isolate(gene_name_list())[label == input$gene,]$value)))
        }
        print(paste("Page:", id, "(General observer)"))
        tx_update_select(isolate(input$gene), gene_name_list, page = id)},
        ignoreNULL = TRUE, ignoreInit = TRUE, priority = -15)
      browser_option_id <-ifelse(id == "browser_allsamp",
                                 "default_gene_meta", "default_gene")

      selected_gene <- ifelse(exists("browser_options"),
                              browser_options[browser_option_id],
                              choices[1])
      gene_update_select_internal(isolate(gene_name_list()), selected_gene,
                                  choices = choices)
    }

    if (uses_libs) {
      observeEvent(libs(), library_update_select(libs),
                   ignoreNULL = TRUE, ignoreInit = FALSE)
      if (id == "browser") {
        observeEvent(input$select_all_btn, {
          print("Pressed select all libs")
          # Ensure that the updateSelectizeInput works properly
          library_update_select(libs, selected = libs())
        }, ignoreNULL = TRUE, ignoreInit = TRUE)  # Ensure the event is triggered on every click, including after initialization
      }
    }
    check_url_for_go_on_init()
    init_round <- FALSE
  }
  )
}

org_and_study_changed_checker <- function(input, output, session) {
  with(rlang::caller_env(), {
    cat("Server startup: "); print(round(Sys.time() - time_before, 2))
    ## Static values
    experiments <- all_exp$name
    ## Set reactive values
    org <- reactiveVal("ALL")
    # Store current and last genome
    exps_dir <- ORFik::config()["exp"]
    df <- reactiveVal(get_exp(browser_options["default_experiment"],
                              experiments, without_readlengths_env, exps_dir))
    df_with <- reactiveVal(get_exp(browser_options["default_experiment"],
                              experiments, with_readlengths_env, exps_dir))
    if (nrow(all_exp_meta) > 0) {
      df_meta <- reactiveVal(get_exp(browser_options["default_experiment_meta"],
                                     all_exp_meta$name, .GlobalEnv, exps_dir))
    } else print("No MegaBrowser exps given, ignoring MegaBrowser exp.")

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
    # gene_name_list <- reactiveVal(names_init)
    gene_name_list <- reactive({
      if(rv$changed == FALSE) {names_init}
      else {get_gene_name_categories(df())}}) %>%
      bindCache(rv$curval) %>%
      bindEvent(rv$changed)
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

allsamples_observer_controller <- function(input, output, session) {
  with(rlang::caller_env(), {
  rv <- reactiveValues(lstval=isolate(df())@txdb,
                       curval=isolate(df())@txdb,
                       genome = "ALL",
                       exp = name(isolate(df())),
                       changed=FALSE)
  observe(if (rv$exp != input$dff & input$dff != "") {
    rv$exp <- input$dff
    message("allsamples browser: dff changed, update rv")
  }) %>%
    bindEvent(input$dff, ignoreInit = TRUE, ignoreNULL = TRUE)

  observe(update_rv_changed(rv), priority = 1) %>%
    bindEvent(rv$curval, ignoreInit = TRUE)
  observe({update_rv(rv, df)}) %>%
    bindEvent(df(), ignoreInit = TRUE)

  observe({df(get_exp(rv$exp, experiments, .GlobalEnv, page = "(allsamples)"))}) %>%
    bindEvent(rv$exp, ignoreInit = TRUE, ignoreNULL = TRUE)

  uses_libs <- FALSE
  org <- reactive("ALL")
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
  study_and_gene_observers(input, output, session)
  })
}

meta_motif_files <- function(df) {
  motif_files <- list.files(file.path(resFolder(df), "meta_collection_tables"), pattern = "fst$", full.names = TRUE)
  names(motif_files) <- sub("\\.fst$", "", gsub(".*_", "", basename(motif_files)))
  return(motif_files)
}

