experiment_update_select <- function(org, all_exp, experiments,
                                     selected = "AUTO") {
  orgs_safe <- if (isolate(org()) == "ALL") {
    unique(all_exp$organism)
  } else isolate(org())
  picks <- experiments[all_exp$organism %in% orgs_safe]
  selected <-
  if (selected == "AUTO") {
    picks[1]
  } else selected
  updateSelectizeInput(
    inputId = "dff",
    choices = picks,
    selected = selected,
    server = FALSE
  )
}

gene_update_select <- function(gene_name_list,
                               selected = gene_name_list()[,2][[1]][1]) {
  updateSelectizeInput(
    inputId = "gene",
    choices = unique(gene_name_list()[,2][[1]]),
    selected = selected,
    server = TRUE
  )
}

gene_update_select_heatmap <- function(gene_name_list, selected = "all") {
  updateSelectizeInput(
    inputId = "gene",
    choices = unique(c(selected, gene_name_list()[,2][[1]])),
    selected = selected,
    server = TRUE
  )
}

tx_update_select <- function(gene = NULL, gene_name_list, additionals = NULL,
                             selected = NULL) {
  print(paste("Updating isoform:"))

  if (is.null(gene)) {
    gene <- isolate(gene_name_list())[value == selected,][1]$label
    if (length(gene) == 0 | is.na(gene)) stop("Isoform does not exist in species!")
  }
  print(paste("Gene set:", gene))
  choices <- gene_name_list()[label == gene,1][[1]]
  choices <- c(additionals, choices)
  if (length(choices) == 0) stop("Gene does not exist in this species")
  if (is.null(selected)) selected <- choices[1]
  if (length(selected) > 1) {
    print(isolate(gene_name_list())[value == selected,][1])
  } else if (selected != "all") print(selected)
  updateSelectizeInput(
    inputId = "tx",
    choices = choices,
    selected = selected,
    server = TRUE
  )
}

frame_type_update_select <- function(selected) {
  updateSelectizeInput(
    inputId = "frames_type",
    label = "Select frames display type",
    choices = c("lines", "columns", "stacks", "area", "heatmap"),
    selected = selected
  )
}

library_update_select <- function(libs) {
  updateSelectizeInput(
    inputId = "library",
    choices = libs(),
    selected = libs()[1],
    server = TRUE
  )
}

condition_update_select <- function(cond) {
  updateSelectizeInput(
    inputId = "condition",
    choices = cond(),
    selected = cond()[1]
  )
}

kmer_update_select <- function(select) {
  updateSliderInput(inputId = "kmer", value = select)
}
