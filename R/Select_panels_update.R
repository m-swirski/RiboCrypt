experiment_update_select <- function(org, all_exp, experiments,
                                     selected = "AUTO") {
  if (!is.null(org)) {
    org <- isolate(org())
  }
  experiment_update_select_isolated(org, all_exp, experiments, selected)
}
experiment_update_select_isolated <- function(org, all_exp, experiments,
                                     selected = "AUTO") {
  if (isTruthy(org)) {
    orgs_safe <- if (isolate(org) == "ALL") {
      unique(all_exp$organism)
    } else isolate(org)
  } else orgs_safe <- unique(all_exp$organism)

  picks <- experiments[all_exp$organism %in% orgs_safe]

  selected <-
    if (!isTruthy(selected) || selected == "AUTO") {
      picks[1]
    } else selected
  updateSelectizeInput(
    inputId = "dff",
    choices = picks,
    selected = selected,
    server = TRUE
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
  cat(paste0("Updating ", id, ": '", selected, "'"), sep = "\n")
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
  tx_update_select_isolated(gene, gene_name_list(), additionals, selected, page)
}

gene_choices_from_gene_list <- function(gene_name_list) {
  unique(gene_name_list[, 2][[1]])
}

gene_exists_in_gene_list <- function(gene_name_list, gene) {
  isTruthy(gene) && gene %in% gene_choices_from_gene_list(gene_name_list)
}

resolve_gene_selection <- function(gene_name_list, preferred = NULL,
                                   fallback = NULL) {
  choices <- gene_choices_from_gene_list(gene_name_list)
  if (length(choices) == 0) return(character())

  candidates <- unique(c(preferred, fallback))
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  selected <- candidates[candidates %in% choices][1]

  if (length(selected) == 0 || is.na(selected)) choices[1] else selected
}

resolve_tx_selection <- function(gene_name_list, gene, preferred = NULL,
                                 fallback = NULL, additionals = NULL) {
  if (!identical(gene, "all") && !gene_exists_in_gene_list(gene_name_list, gene)) {
    return(character())
  }

  isoforms <- tx_from_gene_list(gene_name_list, gene = gene,
                                additionals = additionals)
  if (length(isoforms) == 0) return(character())

  candidates <- unique(c(preferred, fallback))
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  selected <- candidates[candidates %in% isoforms][1]

  if (length(selected) == 0 || is.na(selected)) isoforms[1] else selected
}

browser_default_option_names <- function(id) {
  collection_ids <- c("browser_allsamp", "browser_obs", "selector")
  if (id %in% collection_ids) {
    c(gene = "default_gene_meta", tx = "default_isoform_meta")
  } else {
    c(gene = "default_gene", tx = "default_isoform")
  }
}

url_query_has_value <- function(query, name) {
  value <- query[[name]]
  !is.null(value) && length(value) > 0 && !is.na(value[1]) && nzchar(value[1])
}

gene_for_url_tx <- function(gene_name_list, tx) {
  if (!isTruthy(tx)) return(character())
  tx_match <- gene_name_list[value == tx, label][1]
  if (length(tx_match) == 0 || is.na(tx_match)) character() else tx_match
}

resolve_browser_default_options <- function(browser_options, gene_name_list, id,
                                            query = list()) {
  option_names <- browser_default_option_names(id)
  default_gene <- as.character(browser_options[option_names["gene"]])
  default_tx <- as.character(browser_options[option_names["tx"]])
  query_has_gene <- url_query_has_value(query, "gene")
  query_has_tx <- url_query_has_value(query, "tx")

  preferred_gene <- if (!query_has_gene && query_has_tx) {
    gene_for_url_tx(gene_name_list, default_tx)
  } else {
    default_gene
  }
  selected_gene <- resolve_gene_selection(gene_name_list, preferred = preferred_gene)
  selected_tx <- resolve_tx_selection(
    gene_name_list, selected_gene, preferred = default_tx
  )

  if (!query_has_gene && length(selected_gene) > 0 && isTruthy(selected_gene)) {
    browser_options[option_names["gene"]] <- selected_gene
  }
  if (!query_has_tx && length(selected_tx) > 0 && isTruthy(selected_tx)) {
    browser_options[option_names["tx"]] <- selected_tx
  }
  browser_options
}

tx_update_select_isolated <- function(gene = NULL, gene_name_list, additionals = NULL,
                             selected = NULL, page = "") {
  page <- paste0("(", page, ")")
  isoforms <- tx_from_gene_list(gene_name_list, gene, selected,
                                additionals, page)

  if (is.null(selected)) selected <- isoforms[1]
  if (length(selected) > 1) {
    print(isolate(gene_name_list)[value == selected,][1])
  } else if (selected != "all") print(selected)
  cat(paste0("Updating isoform ", page, ": '", selected, "'"), sep = "\n")
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
                              additionals = NULL, page = "") {

  if (is.null(gene)) {
    gene <- gene_name_list[value == selected,][1]$label
    if (length(gene) == 0 | is.na(gene))
      stop("Isoform does not exist in species!", page)
  } else if (gene == "all") {
    return(c(gene, additionals))
  }
  print(paste("Gene set:", gene))
  isoforms <- gene_name_list[label == gene, 1][[1]]
  isoforms <- c(additionals, isoforms)
  if (length(isoforms) == 0)
    stop("Gene does not exist in this species", page)
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
  library_update_select_safe(libs(), selected, id)
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
  factor_has_2_levels <- length(unique(cond())) > 1
  selected <- ifelse(factor_has_2_levels, 2, 1)

  contrast_levels <- unique(cond())[seq(selected)]
  updateSelectizeInput(
    inputId = "condition",
    choices = cond(),
    selected = contrast_levels
  )
}

kmer_update_select <- function(select) {
  updateSliderInput(inputId = "kmer", value = select)
}
