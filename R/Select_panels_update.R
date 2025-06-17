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
                               selected = choices[1],
                               id = "gene",
                               choices = unique(gene_name_list()[,2][[1]]),
                               server = TRUE) {

  gene_update_select_internal(gene_name_list(),
                              selected = selected,
                              id = id,
                              choices = choices,
                              server = server)
}

gene_update_select_internal <- function(gene_name_list,
                                        selected = choices[1],
                                        id = "gene",
                                        choices = unique(gene_name_list[,2][[1]]),
                                        server = TRUE) {
  print(paste("Updating", paste0(id, ":"), selected))
  updateSelectizeInput(
    inputId = id,
    choices = choices,
    selected = selected,
    server = server
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
                             selected = NULL, page = "") {
  page <- paste0("(", page, ")")
  print(paste("Updating isoform", page, ":"))
  isoforms <- tx_from_gene_list(isolate(gene_name_list()), gene, selected,
                                additionals)

  if (is.null(selected)) selected <- isoforms[1]
  if (length(selected) > 1) {
    print(isolate(gene_name_list())[value == selected,][1])
  } else if (selected != "all") print(selected)
  updateSelectizeInput(
    inputId = "tx",
    choices = isoforms,
    selected = selected,
    server = TRUE
  )
}

motif_update_select <- function(motif_name_list, selected = "") {
  updateSelectizeInput(
    inputId = "motif",
    choices = c(selected, motif_name_list),
    selected = NULL,
    server = TRUE
  )
}

tx_from_gene_list <- function(gene_name_list, gene = NULL, selected = NULL,
                              additionals = NULL) {

  if (is.null(gene)) {
    gene <- gene_name_list[value == selected,][1]$label
    if (length(gene) == 0 | is.na(gene)) stop("Isoform does not exist in species!")
  } else if (gene == "all") {
    return(c(gene, additionals))
  }
  print(paste("Gene set:", gene))
  isoforms <- gene_name_list[label == gene, 1][[1]]
  isoforms <- c(additionals, isoforms)
  if (length(isoforms) == 0) stop("Gene does not exist in this species")
  return(isoforms)
}

frame_type_update_select <- function(selected, id = "frames_type") {
  updateSelectizeInput(
    inputId = id,
    choices = c("lines", "columns", "stacks", "area", "heatmap"),
    selected = selected
  )
}

library_update_select <- function(libs, selected = isolate(libs()[1]),
                                  id = "library") {
  updateSelectizeInput(
    inputId = id,
    choices = libs(),
    selected = selected,
    server = TRUE
  )
}

library_update_select_safe <- function(libs, selected = libs[1],
                                  id = "library") {
  updateSelectizeInput(
    inputId = id,
    choices = libs,
    selected = selected,
    server = TRUE
  )
}

factor_update_select <- function(factor) {
  updateSelectizeInput(
    inputId = "factor",
    choices = factor(),
    selected = factor()[1]
  )
}

condition_update_select <- function(cond) {
  selected <- if (length(unique(cond())) > 1) {
    2
  } else 1
  selected <- unique(cond())[seq(selected)]
  updateSelectizeInput(
    inputId = "condition",
    choices = cond(),
    selected = selected
  )
}

kmer_update_select <- function(select) {
  updateSliderInput(inputId = "kmer", value = select)
}
