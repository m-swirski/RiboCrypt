#' Select box for organism
#'
#' @param genomes name of genomes, returned from list.experiments()
#' @param ns the ID, for shiny session
#' @importFrom Biostrings head
#' @return selectizeInput object
#' @importFrom shinyhelper helper
organism_input_select <- function(genomes, ns) {
  #genomes <- stringr::str_to_sentence(gsub("_", " ", list.genomes()$name))
  selectizeInput(
    inputId = ns("genome"),
    label = "Select an Organism",
    choices = genomes,
    selected = head(genomes, 1),
    multiple = FALSE
  ) %>%
    helper(onclick = "fakeClick('tutorial', 'exp')")
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

gene_input_select <- function(ns, multiple = FALSE) {
    selectizeInput(
    inputId = ns("gene"),
    choices = NULL,
    selected = NULL,
    label = "Select a gene",
    multiple = multiple,
    options = list(placeholder = 'Insert valid gene')
  ) %>%
    helper(onclick = "fakeClick('tutorial', 'gene')")
}

tx_input_select <- function(ns, multiple = FALSE) {
  selectizeInput(
    inputId = ns("tx"),
    choices = "",
    selected = "",
    label = "Select a transcript",
    multiple = multiple,
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
  ) %>%
    helper(onclick = "fakeClick('tutorial', 'lib')")
}

frame_type_select <- function(ns, name = "frames_type",
                              label = "Select frames display type") {
  selectizeInput(
    inputId = ns(name),
    label = label,
    choices = c("lines", "columns", "stacks", "area", "heatmap"),
    selected = "lines",
    multiple = FALSE
  ) %>%
    helper(onclick = "fakeClick('tutorial', 'linetype')")
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

codon_score_input_select <- function(ns) {
  selectizeInput(
    inputId = ns("codon_score"),
    label = "Codon score",
    choices =
      c("percentage", "dispersion(NB)", "alpha(DMN)", "sum"),
    selected = "Percentage",
    multiple = FALSE
  )
}

codon_filter_input_select <- function(ns) {
  numericInput(
    inputId = ns("codon_filter_value"),
    label = "Codon filter value",
    value = 1000,
    min = 0,
    max = NA,
    step = NA
  ) %>%
    helper(onclick = "fakeClick('tutorial', 'codonfilt')")
}

export_format_of_plot <- function(ns) {
  selectizeInput(
    inputId = ns("plot_export_format"),
    label = "Export format",
    choices =
      c("svg", "png"),
    selected = "svg",
    multiple = FALSE
  )
}

condition_input_select <- function(ns, multiple = TRUE) {
  selectizeInput(
    inputId = ns("condition"),
    label = "Select first condition",
    choices = "",
    selected = "",
    multiple = multiple
  ) %>%
    helper(onclick = "fakeClick('tutorial', 'lib')")
}

