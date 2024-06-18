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

observed_cds_heatmap <- function(gene, cds, length_table = NULL,
                                 longestPerGene = TRUE,
                                 minFiveUTR = NULL) {
  if (gene %in% c("", "NULL")) {
    cds()[0]
  } else {
    if ("all" %in% gene) {
      if (!is.null(length_table) & longestPerGene) {
        if (is.null(minFiveUTR)) {
          cds()[length_table()[!duplicated(gene_id) & cds_len > 0,]$tx_name]
        } else {
          cds()[length_table()[!duplicated(gene_id) & cds_len > 0 &
                                 utr5_len >= minFiveUTR,]$tx_name]
        }

      } else cds()

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
  if (mainPlotControls()$region == "Start codon") {
    region <- groupGRangesBy(
      startSites(mainPlotControls()$cds_display, TRUE, TRUE, TRUE))
  } else {
    region <- groupGRangesBy(
      stopSites(mainPlotControls()$cds_display, TRUE, TRUE, TRUE))
  }
  return(region)
}

observed_uorf_annotation <- function(gene, df, all = TRUE, add_uorfs = FALSE) {
  if (add_uorfs) {
    if (all) {
      uorfs <- loadRegion(df(), "uorf")
    } else {
      uorfs <- loadRegion(df(), "uorf", names.keep = gene)
    }
    if (length(uorfs) > 0) names(uorfs) <- paste0("U", seq(length(uorfs)))
    uorfs
  } else GRangesList()
}

observed_translon_annotation <- function(gene, df, all = TRUE, add_translons = FALSE) {
  if (add_translons) {
    translon_annotation <- file.path(dirname(df@fafile),
                                     "predicted_translons",
                                     "predicted_translons_with_sequence_ranges.rds")

    if (file.exists(translon_annotation)) {
      if (all) {
        translons <- readRDS(translon_annotation)
      } else {
        translons <- readRDS(translon_annotation)
        translons <- translons[names(translons) %in% gene]
      }
      if (length(translons) > 0) names(translons) <- paste0("T", seq(length(translons)))
      return(translons)
    }
  }
  return(GRangesList())
}

update_rv <- function(rv, df) {
  print("updating rv txdb")
  # browser()
  if ((rv$lstval == rv$curval) & rv$lstval == "") {
    print("startup updating both")
    rv$lstval <- rv$curval <- df()@txdb
  }
  rv$lstval <- rv$curval; rv$curval <- df()@txdb
  print(isolate(rv$curval))
  print(isolate(rv$lstval))
}

update_rv_changed <- function(rv) {
  if ((rv$lstval != rv$curval)) {
    print("rv_changed")
    if (rv$changed) {
      rv$changed <- FALSE
    } else rv$changed <- TRUE}
  else if (is.null(rv$changed)) rv$changed <- FALSE
}
