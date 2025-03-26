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
