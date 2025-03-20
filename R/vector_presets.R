
normalizations <- function(use_case = "heatmap") {
  norm <-
  if (use_case == "heatmap") {
    c("transcriptNormalized", "zscore", "sum", "log10sum")
  } else if (use_case == "metabrowser") {
    c("maxNormalized", "transcriptNormalized", "zscore","tpm")
  } else stop("Invalid use case option")
  return(norm)
}
