experiment_update_select <- function(org, all_exp, experiments) {
  orgs_safe <- if (isolate(org()) == "ALL") {
    unique(all_exp$organism)
  } else isolate(org())
  picks <- experiments[all_exp$organism %in% orgs_safe]
  updateSelectizeInput(
    inputId = "dff",
    choices = picks,
    selected = picks[1],
    server = FALSE
  )
}

gene_update_select <- function(gene_name_list) {
  updateSelectizeInput(
    inputId = "gene",
    choices = unique(gene_name_list()[,2][[1]]),
    selected = gene_name_list()[,2][[1]][1],
    server = TRUE
  )
}

gene_update_select_heatmap <- function(gene_name_list) {
  updateSelectizeInput(
    inputId = "gene",
    choices = c("all", unique(gene_name_list()[,2][[1]])),
    selected = "all",
    server = TRUE
  )
}

tx_update_select <- function(gene, gene_name_list, additionals = NULL) {
  #browser()
  print("Updating isoform")
  print(isolate(gene_name_list())[label == gene,][1])
  choices <- gene_name_list()[label == gene,1][[1]]
  choices <- c(additionals, choices)
  selected <- choices[1]
  updateSelectizeInput(
    inputId = "tx",
    choices = choices,
    selected = selected,
    server = TRUE
  )
}

library_update_select <- function(libs) {
  updateSelectizeInput(
    inputId = "library",
    choices = libs(),
    selected = libs()[1]
  )
}
