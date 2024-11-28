
#' Load a ORFik collection table
#' @param path the path to gene counts
#' @return a data.table in long format
#' @importFrom fst read_fst
load_collection <- function(path) {
  table  <- fst::read_fst(path)
  setDT(table)
  table[, position := 1:.N, by = library]
  table[, `:=`(library, factor(library, levels = unique(library), ordered = TRUE))]
  return(table)
}

#' Normalize collection table
#'
#' @param table a data.table in long format
#' @param normalization a character string, which mode, for options see RiboCrypt:::normalizations
#' @param lib_sizes named integer vector, default NULL. If given will do a pre tpm normalization
#' for full library sizes
#' @param kmer integer, default 1L (off), if > 1 will smooth out signal with sliding window size kmer.
#' @param add_logscore logical, default TRUE, adds a log(score + 1) to table
#' @param split_by_frame logical, default FALSE
#'  For kmer sliding window, should it split by frame
#' @return a data.table of normalized results
normalize_collection <- function(table, normalization, lib_sizes = NULL,
                                 kmer = 1L, add_logscore = TRUE,
                                 split_by_frame = FALSE) {
  # Sliding window
  if (kmer > 1) table <- smoothenMultiSampCoverage(table, kmer = kmer,
                                                   split_by_frame = split_by_frame)
  # Make tpm
  if (!is.null(lib_sizes)) {
    if (is.character(lib_sizes)) lib_sizes <- readRDS(lib_sizes)
    table[, score_tpm := ((count * 1000)  / lib_sizes[as.integer(library)]) * 10^6]
  } else table[, score_tpm := count]
  # Transcript normalization mode
  norm_opts <- normalizations("metabrowser")
  if (normalization == norm_opts[1]) {
    table[,score := score_tpm / sum(score_tpm), by = library]
    table[,score := score * 1e6]
  } else if (normalization == norm_opts[2]) {
    table[,score := score_tpm / max(score_tpm), by = library]
  } else if (normalization == norm_opts[3]) {
    table[, score := (score_tpm - mean(score_tpm)) / sd(score_tpm), by = library]
  } else table[, score := score_tpm]
  table[is.na(score), score := 0]
  if (add_logscore) table[,logscore := log(score + 1)]
  return(table)
}

match_collection_to_exp <- function(metadata, df) {
  matchings <- chmatch(runIDs(df), metadata$Run)
  matchings <- matchings[!is.na(matchings)]
  if (length(matchings) != nrow(df))
    stop("Metadata does not contain information on all collection samples!")
  return(matchings)
}

#' Cast a collection table to wide format
#' @param table a data.table in long format
#' @param value.var which column to use as scores, default "logscore"
#' @return a table in wide format
collection_to_wide <- function(table, value.var = "logscore") {
  # Remove columns not to be casted to wide format
  # table[, score_tpm := NULL]
  # table[, score := NULL]
  # table[, count := NULL]
  # To wide format
  dtable <- dcast(table, position ~ library, value.var = value.var)
  dtable[, position := NULL]
  return(dtable)
}

#' Get collection table normalized in wide format
#'
#' @inheritParams load_collection
#' @inheritParams normalize_collection
#' @inheritParams collection_to_wide
#' @param df the ORFik experiment to load the precomputed collection from.
#'  It must also have defined runIDs() for all samples.
#' @param metadata a data.table of metadata, must contain the Run column to
#' select libraries.
#' @param metadata_field the column name in metadata, to select to group on.
#' @param as_list logical, default FALSE. Return as list of size 2,
#' count data.table and metadata data.table Set to TRUE if you need metadata
#' subset (needed if you subset the table, to get correct matching)
#' @param subset numeric vector, positional interval to subset, must be <= size of
#' whole region.
#' @param group_on_tx_tpm numeric vector, default NULL.
#' tpm values per libraries. Either for that gene or some other gene.
#' @param min_count integer, default 0. Minimum counts of coverage over transcript
#' to be included.
#' @param format character, default "wide", alternative "long". The format of
#' the table output.
#' @param ratio_interval numeric vector of size 2 or 4, default NULL.
#' If 2, means you should
#' sort libraries on coverage in that region. If 4, means to sort on ratio
#' of that region in this gene vs the other region in another gene.
#' @return a data.table in long or wide (default) format, if as list, it is a
#' list of size 2 (see argument as_list)
compute_collection_table <- function(path, lib_sizes, df,
                                     metadata_field, normalization,
                                     kmer, metadata, min_count = 0, format = "wide",
                                     value.var = "logscore", as_list = FALSE,
                                     subset = NULL, group_on_tx_tpm = NULL,
                                     split_by_frame = FALSE,
                                     ratio_interval = NULL) {
  table <- load_collection(path)
  if (!is.null(subset)) {
    table <- subset_fst_by_interval(table, subset)
  }
  # Normalize
  if (min_count > 0) {
    lib_names <- unique(table$library)
    filt_libs <- table[,.(count = sum(count)),library][count >= min_count,]$library
    table <- table[library %in% filt_libs]
    if (length(filt_libs) == 0)
      stop("Count filter too strict, no libraries with that much reads for this transcript!")
  }
  table <- normalize_collection(table, normalization, lib_sizes, kmer, TRUE,
                                split_by_frame)
  ## # Sort table by metadata column selected
  # Match metadata table and collection runIDs
  if (!is.null(ratio_interval)) {
    tpm <- subset_fst_interval_sum(ratio_interval[seq(2)], table)
    is_ratio <- length(ratio_interval) == 4
    if (is_ratio) {
      tpm2 <- subset_fst_interval_sum(ratio_interval[seq(3,4)], table)
      tpm <- (tpm + 1) / (tpm2 + 1) # Pseudo ratio
    }

    meta_sub <- tpm
    table[, library := factor(library, levels = names(tpm), ordered = TRUE)]
  } else if (!is.null(group_on_tx_tpm)) {
    isoform <- group_on_tx_tpm
    table_path_other <- collection_path_from_exp(df, isoform)
    table_other <- load_collection(table_path_other)
    table_other <- normalize_collection(table_other, "tpm", lib_sizes, 1)
    counts <- table_other[ , .(tpm = sum(score)), by = library]
    tpm <- counts$tpm
    names(tpm) <- counts$library
    meta_sub <- tpm
  } else {
    matchings <- match_collection_to_exp(metadata, df)
    meta_sub <- metadata[matchings, metadata_field, with = FALSE][[1]]
  }
  meta_order <- order(meta_sub)

  table[, library := factor(library, levels = levels(library)[meta_order], ordered = TRUE)]
  if (min_count > 0 & is.null(ratio_interval)) {
    meta_sub <- meta_sub[meta_order][lib_names[meta_order] %in% filt_libs]
  } else if (!is.null(ratio_interval)) {
    meta_sub <- meta_sub[meta_order]
  }
  # Cast to wide format and return
  if (format == "wide") {
    table <- collection_to_wide(table, value.var = value.var)
  }
  if (as_list) return(list(table = table, metadata_field = meta_sub))
  return(table)
}

subset_fst_by_region <- function(df_all, table, id,
                                 gene_mrna = loadRegion(df_all, names.keep = id),
                                 subset = loadRegion(df_all,part = "cds", names.keep = id),
                                 flank_cutoff = 0) {
  stopifnot(is(table, "data.table"))


  subset_pos <- subset_coordinates_grl_to_ir(gene_mrna, subset, flank_cutoff)
  is_long_format <- all(c("library", "position") %in% colnames(table))
  if (is_long_format) {
    return(table[position %in% subset_pos,])
  }
  return(table[subset_pos,])
}

subset_fst_interval_sum <- function(ratio_interval, table) {
  stopifnot(is.numeric(ratio_interval) && length(ratio_interval) == 2)
  stopifnot(all(is.finite(ratio_interval)))

  if (ratio_interval[1] > ratio_interval[2]) stop("Ratio interval must start on >= 1")
  if (ratio_interval[1] < 1) stop("Ratio interval must start on >= 1")
  if (ratio_interval[2] > max(table$position)) stop("Ratio interval must end on <= ncol(heatmap)")

  counts <- table[position %in% seq.int(ratio_interval[1], ratio_interval[2]),
                  .(tpm = sum(score)), by = library]
  tpm <- counts$tpm
  names(tpm) <- counts$library
  return(tpm)
}


subset_fst_coord_by_region <- function(df, id, region_type) {
  extend <- 650 # For yeast
  subset <-
    if (region_type != "mrna") {

      if (organism(df) == "Saccharomyces cerevisiae") {
        region2 <- loadRegion(df, part = "cds", names.keep = id)
        gene_mrna <- extendTrailers(extendLeaders(region2, extend), extend)
      } else gene_mrna <- loadRegion(df, part = "mrna", names.keep = id)

      if (region_type == "leader+cds") {
        if (organism(df) == "Saccharomyces cerevisiae") {
          region <- extendLeaders(GRangesList(startSites(region2, TRUE, FALSE, FALSE)), extend)
        } else {
          region2 <- loadRegion(df, part = "cds", names.keep = id)
          region <- loadRegion(df, part = "leaders", names.keep = id)
        }

        region <- unlistGrl(c(region, region2))
        region <- GRangesList(region)
        names(region) <- id
        subset_coordinates_grl_to_ir(df, id = id, gene_mrna = gene_mrna,
                                     subset = region)
      } else {
        if (organism(df) == "Saccharomyces cerevisiae" & region_type != "cds") {
          if (region_type == "leader") {
            region <- extendLeaders(GRangesList(startSites(region2, TRUE, FALSE, FALSE)), extend)
          } else if (region_type == "trailer") {
            region <- extendTrailers(GRangesList(stopSites(region2, TRUE, FALSE, FALSE)), extend)
          }
        } else region <- loadRegion(df, part = region_type, names.keep = id)
        subset_coordinates_grl_to_ir(df, id = id, gene_mrna = gene_mrna,
                                     subset = region)
      }
    }
  return(subset)
}


subset_coordinates_grl_to_ir <- function(df, id,
                                         gene_mrna = loadRegion(df, names.keep = id),
                                         subset = loadRegion(df,part = "cds", names.keep = id),
                                         flank_cutoff = 0) {
  stopifnot(flank_cutoff >= 0)
  stopifnot(length(subset) == 1)
  stopifnot(length(gene_mrna) == 1)

  subset_txcoord <- pmapToTranscriptF(subset, gene_mrna)
  subset_pos <- seq(as.integer(start(subset_txcoord)) + flank_cutoff, as.integer(end(subset_txcoord)) - flank_cutoff)
  return(subset_pos)
}

subset_fst_by_interval <- function(table, subset) {
  stopifnot(is.numeric(subset) && length(subset) > 0)
  stopifnot(is(table, "data.table"))

  is_long_format <- all(c("library", "position") %in% colnames(table))
  if (is_long_format) {
    return(table[position %in% subset,])
  }
  return(table[subset,])
}

#' Get collection directory
#' @param df ORFik experiment
#' @param must_exists logical, stop if dir does not exists
#' @return file.path(resFolder(df), "collection_tables")
#' @export
collection_dir_from_exp <- function(df, must_exists = FALSE) {
  table_dir <- file.path(resFolder(df), "collection_tables")

  if (must_exists & !file.exists(table_dir))
    stop("There is no collection fst tables directory for this organism,",
         " see vignette for more information on how to make these.")
  return(table_dir)
}

#' Get collection path
#'
#' For directory and id, must be fst format file
#' @param df ORFik experiment
#' @param id character, transcript ids
#' @param gene_name_list a data.table, default NULL, with gene ids
#' @param must_exists logical, stop if dir does not exists
#' @param collection_dir = collection_dir_from_exp(df, must_exists)
#' @return file.path(resFolder(df), "collection_tables")
#' @export
collection_path_from_exp <- function(df, id, gene_name_list = NULL,
                                     must_exists = TRUE,
                                     collection_dir = collection_dir_from_exp(df, must_exists)) {
  table_path <- file.path(collection_dir, paste0(id, ".fst"))
  if (must_exists && !file.exists(table_path)) {
    all_ids_print <- "None"
    if (!is.null(gene_name_list)) {
      all_ids <- gene_name_list[label == gene_name_list[value == id,]$label,]$value
      all_ids_paths <- file.path(collection_dir, paste0(all_ids, ".fst"))
      all_ids <- all_ids[file.exists(all_ids_paths)]
      if (length(all_ids) > 0) {
        all_ids_print <- paste(all_ids, collapse = ", ")
      }
    }
    stop("Gene isoform has no precomputed table, existing isoforms: ", all_ids_print)
  }
  return(table_path)
}
