
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
