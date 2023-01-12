#' Select box for organism
#'
#' @param genomes name of genomes, returned from list.experiments()
#' @param ns the ID, for shiny session
#' @importFrom Biostrings head
#' @return selectizeInput object
organism_input_select <- function(genomes, ns) {
  #genomes <- stringr::str_to_sentence(gsub("_", " ", list.genomes()$name))
  selectizeInput(
    inputId = ns("genome"),
    label = "Select an organism",
    choices = genomes,
    selected = head(genomes, 1),
    multiple = FALSE
  )
}

experiment_input_select <- function(names, ns) {
  selectizeInput(
    inputId = ns("dff"),
    label = "Select an experiment",
    choices = names,
    selected = NULL,
    multiple = FALSE
  )
}

gene_input_select <- function(ns) {
    selectizeInput(
    inputId = ns("gene"),
    choices = NULL,
    selected = NULL,
    label = "Select a gene",
    multiple = FALSE,
    options = list(placeholder = 'Insert valid gene')
  )
}

tx_input_select <- function(ns) {
  selectizeInput(
    inputId = ns("tx"),
    choices = "",
    selected = "",
    label = "Select a transcript",
    multiple = FALSE,
    options = list(placeholder = 'Insert valid tx')
  )
}

library_input_select <- function(ns, multiple = TRUE) {
  selectizeInput(
    inputId = ns("library"),
    label = "Select libraries",
    choices = "",
    selected = "",
    multiple = multiple
  )
}

frame_type_select <- function(ns, name = "frames_type",
                              label = "Select frames display type") {
  selectizeInput(
    inputId = ns(name),
    label = label,
    choices = c("lines", "columns", "stacks", "area", "heatmap"),
    selected = "lines",
    multiple = FALSE
  )
}

heatmap_region_select <- function(ns) {
  selectizeInput(
    inputId = ns("region"),
    label = "View region",
    choices = c("Start codon", "Stop codon"),
    selected = "Start codon",
    multiple = FALSE
  )
}

normalization_input_select <- function(ns) {
  selectizeInput(
    inputId = ns("normalization"),
    label = "Normalization",
    choices =
      c("transcriptNormalized", "zscore", "sum", "log10sum"),
    selected = "transcriptNormalized",
    multiple = FALSE
  )
}
