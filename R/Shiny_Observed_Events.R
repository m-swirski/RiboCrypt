observed_cds <- function(gene, cds) {
  if (gene %in% c("", "NULL")) {
    cds()[0]
  } else {
    if (gene %in% names(cds())) {
      cds()[gene]
    } else cds()[0]
  }
}

observed_cds_heatmap <- function(gene, cds) {
  if (gene %in% c("", "NULL")) {
    cds()[0]
  } else {
    if ("all" %in% gene) {
      cds()
    } else if (gene %in% names(cds())) {
      cds()[gene]
    } else cds()[0]
  }
}


observed_gene <- function(gene, tx) {
  if (gene %in% c("", "NULL")) {
    names(tx()[1])
  } else gene
}

observed_gene_heatmap <- function(gene, tx) {
  if (gene %in% c("", "NULL")) {
    "all"
  } else gene
}

observed_exp_subset <- function(library, libs, df) {
  libs_to_pick <- if (is.null(library)) {
    libs()[1]
  } else library
  df()[which(libs() %in% libs_to_pick),]
}

observed_cds_point <- function(mainPlotControls) {
  class(mainPlotControls()$region)
  if (mainPlotControls()$region == "Start codon") {
    region <- groupGRangesBy(
      startSites(mainPlotControls()$cds_display, TRUE, TRUE, TRUE))
  } else {
    region <- groupGRangesBy(
      stopSites(mainPlotControls()$cds_display, TRUE, TRUE, TRUE))
  }
  return(region)
}
