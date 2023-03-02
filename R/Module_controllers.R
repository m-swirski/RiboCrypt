### NGLVieweR (protein structures) ###
# TODO: Move as much as possible of protein stuff out of page_browser
module_protein <- function(input, output, session) {
  with(rlang::caller_env(), {
# Setup reactive values needed for structure viewer
dynamicVisible <- reactiveVal(FALSE)
selectedRegion <- reactiveVal(NULL)
selectedRegionProfile <- reactive({
  req(selectedRegion())
  result <- cds()[names(cds()) == selectedRegion()] %>%
    getRiboProfile(mainPlotControls()$reads[[1]]) %>%
    (function (x) { x$count[seq.int(1, length(x$count), 3)] })()
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
# NGL viewer widget
protein_structure_dir <- reactive({
  file.path(dirname(df()@fafile), "protein_structure_predictions")
})
region_dir <- reactive({
  file.path(protein_structure_dir(), selectedRegion())
})
pdb_file <- reactive({
  file.path(region_dir(), "ranked_0.pdb")
})
pdb_file_exists <- reactive(pdb_exists(pdb_file))
output$dynamic <- renderNGLVieweR(
  protein_struct_render(pdb_file_exists, selectedRegionProfile, pdb_file))
# Variable UI logic
output$variableUi <- renderUI(
  protein_struct_plot(selectedRegionProfile, dynamicVisible,
                      pdb_file_exists, session))
  }
  )
}

### NGLVieweR (protein structures) ###
# TODO: Move as much as possible of protein stuff out of page_browser
study_and_gene_observers <- function(input, output, session) {
  with(rlang::caller_env(), {
    if (!exists("all_is_gene", mode = "logical")) all_is_gene <- FALSE
    if (!exists("uses_gene", mode = "logical")) uses_gene <- TRUE

    observe(rv$genome <- input$genome) %>%
      bindEvent(input$genome, ignoreInit = TRUE, ignoreNULL = TRUE)
    observe(rv$exp <- input$dff) %>%
      bindEvent(input$dff, ignoreInit = TRUE, ignoreNULL = TRUE)
    observe(updateSelectizeInput(
      inputId = "genome",
      choices = c("ALL", unique(all_exp$organism)),
      selected = rv$genome,
      server = TRUE
    )) %>%
      bindEvent(rv$genome, ignoreInit = TRUE, ignoreNULL = TRUE)
    observeEvent(rv$exp, experiment_update_select(org, all_exp, experiments,
                                                  rv$exp),
                 ignoreInit = TRUE)

    observeEvent(org(), experiment_update_select(org, all_exp, experiments),
                 ignoreInit = TRUE)
    if (all_is_gene) {
      updateSelectizeInput(
        inputId = "gene",
        choices = c("all", unique(isolate(gene_name_list())[,2][[1]])),
        selected = "all",
        server = TRUE
      )
      observeEvent(gene_name_list(), gene_update_select_heatmap(gene_name_list),
                   ignoreInit = TRUE)
      observeEvent(input$gene, tx_update_select(isolate(input$gene),
                                                gene_name_list, "all"),
                   ignoreNULL = TRUE, ignoreInit = TRUE)
    } else if (uses_gene) {
      # TODO: decide if updateSelectizeInput should be on top here or not
      observeEvent(gene_name_list(), gene_update_select(gene_name_list),
                   ignoreNULL = TRUE, ignoreInit = TRUE)
      observeEvent(input$gene, tx_update_select(isolate(input$gene),
                                                gene_name_list),
                   ignoreNULL = TRUE, ignoreInit = TRUE)
      updateSelectizeInput(
        inputId = "gene",
        choices = unique(isolate(gene_name_list())[,2][[1]]),
        selected = isolate(input$gene),
        server = TRUE
      )
    }
    observeEvent(libs(), library_update_select(libs))
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

    observe(org(rv$genome)) %>%
      bindEvent(rv$genome, ignoreInit = TRUE, ignoreNULL = TRUE)
    observe({df(get_exp(rv$exp, experiments, without_readlengths_env))}) %>%
      bindEvent(rv$exp, ignoreInit = TRUE, ignoreNULL = TRUE)
    observe({df_with(get_exp(rv$exp, experiments, with_readlengths_env))}) %>%
      bindEvent(rv$exp, ignoreInit = TRUE, ignoreNULL = TRUE)

    cat("Pre modules: "); print(round(Sys.time() - time_before, 2))
  }
  )
}



