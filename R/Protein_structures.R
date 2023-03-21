# Utilities or computing coloring scheme based on Ribo-Seq results
valuesToColors <- function(vals) {
  palette <- colorRampPalette(c("blue", "green", "yellow", "orange", "red"), bias = 0.5)(1001)
  palette[(vals / max(vals) * 1000) + 1]
}

protein_struct_render <- function(selectedRegionProfile, structure_file) {
  req(structure_file(), selectedRegionProfile())
  structure_file() %>% NGLVieweR() %>%
    stageParameters(backgroundColor = "white") %>%
    onRender(fetchJS("sequence_viewer_coloring.js"),
             valuesToColors(selectedRegionProfile()))
}

pdb_exists <- function(pdb_file) {
  req(pdb_file())
  file.exists(pdb_file())
}

protein_struct_plot <- function(selectedRegion, selectedRegionProfile, dynamicVisible,
                                session, structureChoices = list()) {
  req(dynamicVisible(), selectedRegionProfile())
  
  ns <- session$ns
  
  
  widgetCloseBtn <- actionButton(ns("dynamicClose"), "Close")
  widgetHeader <- renderText(selectedRegion())
  widgetSelector <- selectInput(ns("structureViewerSelector"), "Select structure", structureChoices())
  
  list(widgetCloseBtn, widgetHeader, widgetSelector, NGLVieweROutput(ns("dynamic")))
}
