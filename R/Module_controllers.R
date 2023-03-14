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
    # NGL viewer widget
    protein_structure_dir <- reactive({
      file.path(dirname(df()@fafile), "protein_structure_predictions")
    })
    region_dir <- reactive({
      file.path(protein_structure_dir(), selectedRegion())
    })
    pdb_files <- reactive(list.files(region_dir()))
    pdb_file <- reactive({
      if(is.null(input$structureViewerSelector)) {
        file.path(region_dir(), head(pdb_files()))
      } else {
        file.path(region_dir(), input$structureViewerSelector)
      }
    })
    pdb_file_exists <- reactive(pdb_exists(pdb_file))
    output$dynamic <- renderNGLVieweR(protein_struct_render(pdb_file_exists, selectedRegionProfile, pdb_file))
    # Variable UI logic
    output$variableUi <- renderUI(
      protein_struct_plot(
        selectedRegion,
        selectedRegionProfile,
        dynamicVisible,
        pdb_file_exists,
        session,
        pdb_files()
      )
    )
  })
}

### NGLVieweR (protein structures) ###
# TODO: Move as much as possible of protein stuff out of page_browser
study_and_gene_observers <- function(input, output, session) {
  with(rlang::caller_env(), {
    if (!exists("all_is_gene", mode = "logical")) all_is_gene <- FALSE
    observeEvent(org(), experiment_update_select(org, all_exp, experiments))
    if (all_is_gene) {
      observeEvent(gene_name_list(), gene_update_select_heatmap(gene_name_list))
      observeEvent(input$gene, tx_update_select(isolate(input$gene),
                                                gene_name_list, "all"),
                   ignoreNULL = TRUE, ignoreInit = TRUE)
    } else {
      observeEvent(gene_name_list(), gene_update_select(gene_name_list))
      observeEvent(input$gene, tx_update_select(isolate(input$gene),
                                                gene_name_list),
                   ignoreNULL = TRUE, ignoreInit = TRUE)
    }
    observeEvent(libs(), library_update_select(libs))
  }
  )
}

org_and_study_changed_checker <- function(input, output, session) {
  with(rlang::caller_env(), {
    # Static values
    experiments <- all_exp$name
    # Set reactive values
    org <- reactive(input$genome)
    rv <- reactiveValues(lstval="",curval="") # Store current and last genome
    rv_changed <- reactiveVal(NULL) # Did genome change?
    observe(update_rv_changed(rv, rv_changed), priority = 1) %>%
      bindEvent(rv$curval)
    df <- reactive(get_exp(input$dff, experiments, env))
    observeEvent(df(), update_rv(rv, df), priority = 2)
    libs <- reactive(bamVarName(df()))
  }
  )
}

