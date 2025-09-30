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
