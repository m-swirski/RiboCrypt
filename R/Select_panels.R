experiment_input_select <- function(names) {
  selectizeInput(
    inputId = "dff",
    label = "Select an experiment",
    choices = names,
    selected = head(names),
    multiple = FALSE
  )
}

gene_input_select <- function() {
    selectizeInput(
    inputId = "gene",
    choices = "",
    selected = "",
    label = "Select a gene",
    multiple = FALSE
  )
}

library_input_select <- function() {
  selectizeInput(
    inputId = "library",
    label = "Select libraries",
    choices = "",
    selected = "",
    multiple = TRUE
  )
}

frame_type_select <- function() {
  selectizeInput(
    inputId = "frames_type",
    label = "Select frames display type",
    choices = c("lines", "columns", "stacks", "area"),
    selected = "lines",
    multiple = FALSE
  )
}
