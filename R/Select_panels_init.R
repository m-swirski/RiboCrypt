#' Select box for organism
#'
#' @param genomes name of genomes, returned from list.experiments()
#' @param ns the ID, for shiny session
#' @importFrom Biostrings head
#' @return selectizeInput object
#' @importFrom shinyhelper helper
organism_input_select <- function(genomes, ns) {
  selectizeInput(
    inputId = ns("genome"),
    label = "Select an Organism",
    choices = genomes,
    selected = head(genomes, 1),
    multiple = FALSE
  ) %>%
    helper(onclick = "fakeClick('tutorial', 'exp')")
}

experiment_input_select <- function(names, ns, browser_options = NULL,
                                    option_name = "default_experiment") {
  selectizeInput(
    inputId = ns("dff"),
    label = "Select an experiment",
    choices = names,
    selected = browser_options[option_name],
    multiple = FALSE
  )
}

gene_input_select <- function(ns, multiple = FALSE, browser_options = NULL,
                              choices = as.character(browser_options["default_gene"]),
                              label = "Select a gene",
                              id = "gene",
                              selected = as.character(browser_options["default_gene"])) {
    selectizeInput(
    inputId = ns(id),
    choices = choices,
    selected = selected,
    label = label,
    multiple = multiple,
    options = list(placeholder = 'Insert valid gene')
  ) %>%
    helper(onclick = "fakeClick('tutorial', 'gene')")
}

tx_input_select <- function(ns, multiple = FALSE, choices = NULL,
                            selected = choices$value[1]) {
  selectizeInput(
    inputId = ns("tx"),
    choices = choices$value,
    selected = selected,
    label = "Select a transcript",
    multiple = multiple,
    options = list(placeholder = 'Insert valid tx')
  )
}

library_input_select <- function(ns, multiple = TRUE, choices = "",
                                 selected = if (is(choices, "data.frame")) {
                                   choices[1, , drop = FALSE]
                                 } else {
                                   choices[1]
                                   },
                                 label = "Select libraries", id = "library") {
  selectizeInput(
    inputId = ns(id),
    label = label,
    choices = choices,
    selected = selected,
    multiple = multiple,
    width = "100%",
    options = list(plugins = list("remove_button"))
  ) %>%
    helper(onclick = "fakeClick('tutorial', 'lib')")
}

motif_input_select <- function(ns, multiple = FALSE, browser_options = NULL,
                              choices = NULL,
                              label = "or Select a meta motif",
                              id = "motif",
                              selected = NULL) {
  selectizeInput(
    inputId = ns(id),
    choices = choices,
    selected = selected,
    label = label,
    multiple = multiple,
    options = list(placeholder = 'Select valid Motif')
  )
}

frame_type_select <- function(ns, name = "frames_type",
                              label = "Select frames display type",
                              selected = "lines") {
  selectizeInput(
    inputId = ns(name),
    label = label,
    choices = c("lines", "columns", "stacks", "area", "heatmap", "animate"),
    selected = selected,
    multiple = FALSE
  ) %>%
    helper(onclick = "fakeClick('tutorial', 'linetype')")
}

frame_subsetter_select <- function(ns, name = "frames_subset",
                                   label = "Subset frames",
                                   selected = NULL) {
  selectizeInput(
    inputId = ns(name),
    label = label,
    choices = c("all", "red", "green", "blue"),
    selected = selected,
    multiple = TRUE
  ) %>%
    helper(onclick = "fakeClick('tutorial', 'settings')")
}

region_view_select <- function(ns, name, label,
                               selected = "mrna") {
  selectizeInput(
    inputId = ns(name),
    label = label,
    choices = c("mrna", "cds", "leader+cds", "leader", "trailer"),
    selected = selected,
    multiple = FALSE
  )
}

heatmap_region_select <- function(ns, add = NULL) {
  choices <- c("Start codon", "Stop codon", add)
  selectizeInput(
    inputId = ns("region"),
    label = "Anchor region",
    choices = choices,
    selected = choices[1],
    multiple = FALSE
  )
}

heatmap_color_select <- function(ns, name = "heatmap_color",
                              label = "Color theme",
                              selected = "Matrix (black,green,red)") {
  selectizeInput(
    inputId = ns(name),
    label = label,
    choices = c("default (White-Blue)", "Matrix (black,green,red)"),
    selected = selected,
    multiple = FALSE
  ) %>%
    helper(onclick = "fakeClick('tutorial', 'mbrowser')")
}

normalization_input_select <- function(ns,
                                       choices = normalizations(),
                                       selected = choices[1],
                                       help_link = "heatmap") {
  selectizeInput(
    inputId = ns("normalization"),
    label = "Normalization",
    choices = choices,
    selected = selected,
    multiple = FALSE
  ) %>%
    helper(onclick = paste("fakeClick('tutorial',", paste0("'", help_link, "')")))
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

codon_filter_input_select <- function(ns, value = 1000) {
  numericInput(
    inputId = ns("codon_filter_value"),
    label = "Codon filter value",
    value = value,
    min = 0,
    max = NA,
    step = NA
  ) %>%
    helper(onclick = "fakeClick('tutorial', 'codonfilt')")
}

diff_method_input_select <- function(ns) {
  selectizeInput(
    inputId = ns("diff_method"),
    label = "Differential method",
    choices =
      c("FPKM ratio", "DESeq2"),
    multiple = FALSE
  ) %>%
    helper(onclick = "fakeClick('tutorial', 'diffexp')")
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
    label = "Select two conditions",
    choices = "",
    selected = "",
    multiple = multiple
  ) %>%
    helper(onclick = "fakeClick('tutorial', 'diffexp')")
}

metadata_input_select <- function(ns, metadata,
                                  multiple = FALSE, browser_options = NULL,
                                  choices = colnames(metadata), selected = "TISSUE",
                                  label = "Order on:", id = "metadata",
                                  add = NULL) {
  choices <- c(choices, add)
  selectizeInput(
    inputId = ns(id),
    choices = choices,
    selected = selected,
    label = label,
    multiple = multiple,
    options = list(placeholder = label)
  )
  #  %>%
  #   helper(onclick = "fakeClick('tutorial', 'metadata')")
}

umap_color_by_input_select <- function(ns, names = "auto",
                                       browser_options = NULL,
                                       option_name = "default_experiment") {
  if (names == "auto") {
    names <- c('Tissue' = "tissue", 'Cell line' = "cell_line",
               'Inhibitor' = "inhibitors", 'BioProject' = "BioProject",
               'Author' = "author")
  }
  selectizeInput(
    inputId = ns("umap_col"),
    label = "Color on",
    choices = names,
    selected = names[seq(2)],
    multiple = TRUE
  )
}

umap_plot_type <- function(ns, names = c("UMAP", "UMAP centroids"),
                           browser_options = NULL,
                           option_name = "default_experiment") {
  selectizeInput(
    inputId = ns("umap_plot_type"),
    label = "Plot type",
    choices = names,
    selected = names[1],
    multiple = FALSE
  )
}

