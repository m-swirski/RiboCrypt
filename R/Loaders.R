load_reads <- function(dff, read_type, BPPARAM = BiocParallel::SerialParam()) {
  force(
    outputLibs(
      dff,
      type = read_type,
      output.mode = "envirlist",
      naming = "full",
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
