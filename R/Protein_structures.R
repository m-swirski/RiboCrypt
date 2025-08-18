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

coordinates_to_pep_id_path <- function(grl_coordinates, linker_dt, pep_dir) {
  stopifnot(is(linker_dt, "data.table"))
  stopifnot(colnames(linker_dt) == c("ID", "coordinates"))
  selected_linker <- linker_dt[coordinates == as.character(grl_coordinates),]
  if (nrow(selected_linker) == 0) {
    print("Translon does not have a peptide linker ID, is there a longer one?")
    return("")
  }
  paths <- pep_id_to_path(selected_linker[1,]$ID, pep_dir)
  return(paths)
}

pep_id_to_path <- function(id, pep_dir) {
  regex_str <- paste0("^", id, "_unrelaxed_rank_001_alphafold2_ptm_.*\\.pdb$")
  pep_dirs <- file.path(pep_dir, c(seq(0:17), paste0("nRNA_", seq(9))), "predictions")
  pep_dirs <- pep_dirs[dir.exists(pep_dirs)]
  paths <- unlist(lapply(pep_dirs, function(x)
    list.files(x, regex_str, full.names = TRUE)), use.names = FALSE)
  if (length(paths) == 0) {
    print("Linker peptide id found, but file does not exist!")
    return("")
  }
  return(paths)
}

protein_struct_plot <- function(selectedRegion, selectedRegionProfile, dynamicVisible,
                                session, structureChoices = list()) {
  req(dynamicVisible(), selectedRegionProfile())
  ns <- session$ns

  widgetCloseBtn <- actionButton(
    ns("dynamicClose"),
    label = icon("times", class = "text-white"),
    class = "btn btn-sm btn-primary",
    style = "width: 100%; background-color: #007bff; border-color: #007bff;",
    title = "Remove protein structure"
  )

  widgetHeader <- tags$div(
    h4(paste("Protein isoform:", selectedRegion())),
    style = "display: flex; align-items: flex-end; height: 100%;"
  )

  widgetSelector <- selectInput(ns("structureViewerSelector"), NULL, structureChoices(), width = "100%")

  tagList(
    tags$div(
      style = "border-top: 2px solid black; margin-top: 15px; padding-top: 10px;",
      fluidRow(
        column(5),  # Blank space for padding left
        column(4, widgetHeader),
        column(2, widgetSelector),
        column(1, widgetCloseBtn)
      )
    ),
    fluidRow(NGLVieweROutput(ns("dynamic")))
  )
}
