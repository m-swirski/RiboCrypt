# Utilities or computing coloring scheme based on Ribo-Seq results
valuesToColors <- function(vals) {
  palette <- colorRampPalette(c("blue", "green", "yellow", "orange", "red"), bias = 0.5)(1001)
  palette[(vals / max(vals) * 1000) + 1]
}

protein_struct_render <- function(pdb_file_exists, selectedRegionProfile,
                                  pdb_file) {
  req(pdb_file_exists(), selectedRegionProfile())
  pdb_file() %>% NGLVieweR() %>%
    stageParameters(backgroundColor = "white") %>%
    onRender(fetchJS("sequence_viewer_coloring.js"),
             valuesToColors(selectedRegionProfile()))
}

pdb_exists <- function(pdb_file) {
  req(pdb_file())
  file.exists(pdb_file())
}

protein_struct_plot <- function(selectedRegionProfile, dynamicVisible,
                                pdb_file_exists, session, structureChoices = list()) {
  req(dynamicVisible(), pdb_file_exists(), selectedRegionProfile())
  
  ns <- session$ns
  
  
  closeBtn <- actionButton(ns("dynamicClose"), "Close")
  structureSelector <- selectInput("structureSelector", "Select structure", structureChoices)
  
  list(closeBtn, structureSelector, NGLVieweROutput(ns("dynamic")))
}
