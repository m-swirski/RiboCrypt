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
  
  
  widgetCloseBtn <- actionButton(ns("dynamicClose"), "Close", width = "100%")
  widgetHeader <- h3(selectedRegion(), style = "text-align: center; transform: translate(0%, -60%);")
  widgetSelector <- selectInput(ns("structureViewerSelector"), NULL, structureChoices(), width = "100%")
  
  tagList(
    fluidRow(
      column(2, widgetCloseBtn),
      column(6, widgetHeader, offset = 1),
      column(2, widgetSelector, offset = 1)
      ),
    fluidRow(NGLVieweROutput(ns("dynamic")))
    )
}
