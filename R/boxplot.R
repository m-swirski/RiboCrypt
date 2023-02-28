distribution_plot <- function(df, display_region = NULL, annotation = NULL, leaders_extension = NULL, trailers_extension = NULL) {
  
  sct <- countTable(df, type = "summarized")
  
  ct <- countTable(df, type = "fpkm")
  names <- rownames(sct)
  
  mvars <- colnames(ct)
  ct$names <- names
  ct <- melt(ct, measure.vars = mvars)
  colnames(ct) <- c("names", "library", "RPKM")
  display_region <- extendTrailers(extendLeaders(display_region, leaders_extension), trailers_extension)
  overlaps <- subsetByOverlaps(annotation, display_region)
  p <- ggplot(ct, aes(x = library, group = library, y  = log2(RPKM))) +
    geom_violin() +
    geom_boxplot(width=0.1) +
    geom_jitter(data = ct[names %in% names(overlaps)], mapping = aes(x = library, fill = names, color = names, y  = log2(RPKM),text = names)) +
    theme_bw()
  return(ggplotly(p))
  
}


# get_gene_name_categories(df)$value