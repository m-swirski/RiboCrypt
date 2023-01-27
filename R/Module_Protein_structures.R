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
