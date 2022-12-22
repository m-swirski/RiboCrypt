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
    selected = "",
    multiple = FALSE
  )
}

gene_input_select <- function(ns) {
    selectizeInput(
    inputId = ns("gene"),
    choices = "",
    selected = "",
    label = "Select a gene",
    multiple = FALSE
  )
}

library_input_select <- function(ns) {
  selectizeInput(
    inputId = ns("library"),
    label = "Select libraries",
    choices = "",
    selected = "",
    multiple = TRUE
  )
}

frame_type_select <- function(ns) {
  selectizeInput(
    inputId = ns("frames_type"),
    label = "Select frames display type",
    choices = c("lines", "columns", "stacks", "area"),
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
