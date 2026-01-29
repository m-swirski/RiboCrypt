get_track_paths <- function(dff) {
  libtypes <- ORFik:::libraryTypes(dff, FALSE)
  must_check_non_unique_mappers_for_non_ribo <- uniqueMappers(dff) &
    !identical("RFP", unique(libtypes))

  if (must_check_non_unique_mappers_for_non_ribo) {
    index_rfp <- which(libtypes == "RFP")
    index_not_rfp <- which(libtypes != "RFP")
    reads <- list()
    if (length(index_rfp) > 0) reads[index_rfp] <- get_track_paths_internal(dff[index_rfp,])
    if (length(index_not_rfp) > 0) {
      dff_non_unique <- dff[index_not_rfp,]
      uniqueMappers(dff_non_unique) <- FALSE
      reads[index_not_rfp] <- get_track_paths_internal(dff_non_unique)
    }
    stopifnot(all(lengths(reads) > 0))
  } else {
    reads <- get_track_paths_internal(dff)
  }

  return(reads)
}

get_track_paths_internal <- function(dff) {
  reads <- try(filepath(dff, "bigwig", suffix_stem = c("_pshifted", "")), silent = TRUE)
  invalid_reads <- is(reads, "try-error") ||
    (!all(file.exists(unlist(reads, use.names = FALSE))) |
       any(duplicated(unlist(reads, use.names = FALSE))))
  if (invalid_reads) {
    reads <- filepath(dff, "bigwig", suffix_stem = c("_pshifted", ""),
                      base_folders = libFolder(dff, "all"))
  }
  return(reads)
}

load_reads <- function(dff, prefered_read_type, validate_libs = FALSE,
                       BPPARAM = BiocParallel::SerialParam()) {
  pref_dir <- if (prefered_read_type == "cov") {
    "cov_RLE"
  } else if (prefered_read_type == "covl") {
    "cov_RLE_List"
  } else stop("Only cov and covRLE supported safely at the moment!")
  read_type <- ifelse(dir.exists(file.path(libFolder(dff), pref_dir)), prefered_read_type,
                      "pshifted")
  message("Using read type: ", read_type)
  paths <- filepath(dff, read_type, suffix_stem = c("_pshifted", ""))
  if (length(paths) > 0) {
    message("First file to load is:")
    paths[1]
  }

  force(
    outputLibs(
      dff,
      type = read_type,
      paths = paths,
      output.mode = "envirlist",
      naming = "fullexp",
      validate_libs = validate_libs,
      BPPARAM = BPPARAM
    )
  )
}

load_custom_regions <- function(useCustomRegions, df) {
  if(isTRUE(useCustomRegions)) {
    protein_structure_path <- file.path(dirname(df()@fafile), "protein_structure_predictions", "custom_regions.csv")
    if (file.exists(protein_structure_path)) {
      orfs_flt <- fread(protein_structure_path)
      orfs_flt_grl <- GRanges(orfs_flt) %>% groupGRangesBy(.,.$names)
    } else NULL
  } else NULL
}

# Function to load data
load_data <- function(species) {

  data <- load_data_internal(species)
  reactiveValues(translon_table = data$translon_table, df = data$df)
}

load_data_internal <- function(species) {
  df <- read.experiment(species, validate = FALSE)
  table_path <- file.path(refFolder(df),
                          "predicted_translons",
                          "predicted_translons_with_sequence.fst")
  if (file.exists(table_path)) {
    translon_table <- fst::read_fst(table_path, as.data.table = TRUE)
    setattr(translon_table, "exp", species)
  } else {
    NULL
  }
  return(list(translon_table = translon_table, df = df))
}

# Function to load data
load_data_umap <- function(species, color.by = NULL) {

  data <- load_data_umap_internal(species, color.by)
  reactiveValues(dt_umap = data$dt_umap, df = data$df)
}

load_data_umap_internal <- function(species, color.by = c("tissue", "cell_line")) {
  df <- read.experiment(species, validate = FALSE)
  dir <- file.path(refFolder(df), "UMAP")
  table_path <- file.path(dir, "UMAP_by_gene_counts.fst")
  if (file.exists(table_path)) {
    dt_umap <- fst::read_fst(table_path, as.data.table = TRUE)
    if (length(color.by) > 1) {
      dt_umap[, color_column := do.call(paste, c(.SD, sep = " | ")), .SDcols = color.by]
    } else dt_umap[, color_column := get(color.by)]

    setattr(dt_umap, "exp", species)
    setattr(dt_umap, "color.by", color.by)
  } else {
    stop("Species has no computed UMAP, pick another!")
  }
  return(list(dt_umap = dt_umap, df = df))
}
