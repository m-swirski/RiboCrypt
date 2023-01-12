observed_cds_annotation <- function(gene, cds, all = TRUE) {
  if (gene %in% c("", "NULL")) {
    cds()[0]
  } else {
    if (all) {
      cds()
    } else {
      if (gene %in% names(cds())) {
        cds()[gene]
      } else cds()[0]
    }

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


observed_tx_annotation <- function(gene, tx) {
  if (gene %in% c("", "NULL") | !(gene %in% names(tx()))) {
    print("Empty tx input, using first tx in organism")
    tx()[1]
  } else tx()[gene]
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

update_rv <- function(rv, df) {
  print("updating rv")
  if ((rv$lstval == rv$curval) & rv$lstval == "") {
    print("startup updating both")
    rv$lstval <- rv$curval <- df()@txdb
  }
  rv$lstval <- rv$curval; rv$curval <- df()@txdb
  print(isolate(rv$curval))
  print(isolate(rv$lstval))
}

update_rv_changed <- function(rv, rv_changed) {
  {if ((rv$lstval != rv$curval)) {print("rv_changed"); if (rv_changed()) {rv_changed(F)} else (rv_changed(T))}
    else if (is.null(rv_changed())) {rv_changed(F) }}
}
