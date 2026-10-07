get_track_paths <- function(dff) {
  libtypes <- ORFik::libraryTypes(dff, FALSE)
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
  get_library_paths_internal(dff, "bigwig")
}

#' Check that one distinct set of files was found for every library.
#' @noRd
library_paths_are_valid <- function(paths, expected_libraries) {
  if (inherits(paths, "try-error") || length(paths) != expected_libraries) {
    return(FALSE)
  }
  flat_paths <- unlist(paths, use.names = FALSE)
  length(flat_paths) >= expected_libraries &&
    all(nzchar(flat_paths)) &&
    all(file.exists(flat_paths)) &&
    !anyDuplicated(flat_paths)
}

#' Resolve derived library files across ordinary and collection experiments.
#' @noRd
get_library_paths_internal <- function(dff, read_type,
                                       suffix_stem = c("_pshifted", "")) {
  paths <- try(
    filepath(dff, read_type, suffix_stem = suffix_stem),
    silent = TRUE
  )
  if (library_paths_are_valid(paths, nrow(dff))) return(paths)

  # ORFik::filepath() normally uses the first library folder for every row.
  # Collections can contain rows from several study folders, so retry with
  # the folder belonging to each individual library.
  paths <- filepath(
    dff,
    read_type,
    suffix_stem = suffix_stem,
    base_folders = libFolder(dff, "all")
  )
  if (!library_paths_are_valid(paths, nrow(dff))) {
    stop("Could not resolve distinct ", read_type,
         " files for every selected library.", call. = FALSE)
  }
  paths
}

#' Resolve covRle files, including experiments spanning multiple folders.
#' @noRd
get_covRle_paths_internal <- function(dff) {
  get_library_paths_internal(dff, "cov")
}

#' Resolve covRleList files, including experiments spanning multiple folders.
#' @noRd
get_covRleList_paths_internal <- function(dff) {
  get_library_paths_internal(dff, "covl")
}

#' Resolve p-shifted files, including fallback from unavailable coverage files.
#' @noRd
get_pshifted_paths_internal <- function(dff) {
  get_library_paths_internal(dff, "pshifted")
}

load_reads <- function(dff, prefered_read_type, validate_libs = FALSE,
                       BPPARAM = BiocParallel::SerialParam()) {
  preferred_path_loader <- switch(
    prefered_read_type,
    cov = get_covRle_paths_internal,
    covl = get_covRleList_paths_internal,
    stop("Only covRle ('cov') and covRleList ('covl') are supported.")
  )
  paths <- try(preferred_path_loader(dff), silent = TRUE)
  preferred_available <- !inherits(paths, "try-error")
  read_type <- if (preferred_available) prefered_read_type else "pshifted"
  if (!preferred_available) paths <- get_pshifted_paths_internal(dff)

  message("Using read type: ", read_type)
  if (length(paths) > 0) {
    message("First file to load is:")
    message(unlist(paths, use.names = FALSE)[[1]])
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

#' Load covRle coverage with collection-aware path resolution.
#' @noRd
load_covRle <- function(dff, validate_libs = FALSE,
                        BPPARAM = BiocParallel::SerialParam()) {
  load_reads(dff, "cov", validate_libs = validate_libs, BPPARAM = BPPARAM)
}

#' Load covRleList coverage with collection-aware path resolution.
#' @noRd
load_covRleList <- function(dff, validate_libs = FALSE,
                            BPPARAM = BiocParallel::SerialParam()) {
  load_reads(dff, "covl", validate_libs = validate_libs, BPPARAM = BPPARAM)
}

#' Load the shift table belonging to one selected collection library.
#' @noRd
load_library_shift_table <- function(dff) {
  if (nrow(dff) != 1L) {
    stop("A single selected library is required to load its shift table.")
  }

  shift_path <- file.path(
    libFolder(dff, "all"), "pshifted", "shifting_table.rds"
  )
  if (!file.exists(shift_path)) {
    warning("Shift table not found for the selected library.")
    return(data.table::data.table())
  }

  shifts <- shifts_load(dff, path = shift_path)
  if (!is.list(shifts) || !length(shifts)) {
    warning("The selected library's shift table is empty or malformed.")
    return(data.table::data.table())
  }

  shift_names <- names(shifts)
  selected_stem <- orfik_remove_file_ext(basename(dff$filepath[[1]]))
  shift_stems <- orfik_remove_file_ext(basename(shift_names))
  matching_shift <- which(shift_stems == selected_stem)
  unnamed_single <- length(shifts) == 1L &&
    (is.null(shift_names) || is.na(shift_names[[1]]) || !nzchar(shift_names[[1]]))
  if (!length(matching_shift) && unnamed_single) matching_shift <- 1L
  if (length(matching_shift) != 1L) {
    warning("Could not match the selected library to its study shift table.")
    return(data.table::data.table())
  }

  shift_table <- data.table::as.data.table(shifts[[matching_shift]])
  required_columns <- c("fraction", "offsets_start")
  if (!all(required_columns %in% names(shift_table))) {
    warning("The selected library's shift table is malformed.")
    return(data.table::data.table())
  }
  shift_table
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
  df <- rc_read_experiment(species, validate = FALSE)
  rc_access_reference(df)
  translon_table <- NULL
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
}

load_data_umap_internal <- function(species, color.by = c("tissue", "cell_line")) {
  df <- rc_read_experiment(species, validate = FALSE)
  dir <- file.path(refFolder(df), "UMAP")
  table_path <- file.path(dir, "UMAP_by_gene_counts.fst")
  if (file.exists(table_path)) {
    dt_umap <- fst::read_fst(table_path, as.data.table = TRUE)
    dt_umap <- rc_access_umap(dt_umap, df)
    if (length(color.by) > 1) {
      dt_umap[, color_column := do.call(paste, c(.SD, sep = " | ")), .SDcols = color.by]
    } else dt_umap[, color_column := get(color.by)]

    setattr(dt_umap, "exp", species)
    setattr(dt_umap, "color.by", color.by)
  } else {
    stop("Species has no computed UMAP, pick another!")
  }
  dt_umap
}
