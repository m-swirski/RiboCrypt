distribution_plot <- function(df, display_region = NULL, annotation = NULL, leaders_extension = NULL, trailers_extension = NULL) {
  
  sct <- countTable(df, type = "summarized")
  
  ct <- countTable(df, type = "fpkm")
  names <- rownames(sct)
  
  mvars <- colnames(ct)
  ct$transcript <- names
  ct <- melt(ct, measure.vars = mvars)
  colnames(ct) <- c("transcript", "library", "RPKM")
  display_region <- extendTrailers(extendLeaders(display_region, leaders_extension), trailers_extension)
  display_range <- display_region
  #boilerplate from createGeneModelPanel
  overlaps <- subsetByOverlaps(annotation, display_range)
  plot_width <- widthPerGroup(display_range)
  onames <- rep(names(overlaps), numExonsPerGroup(overlaps, FALSE))
  overlaps <- unlistGrl(overlaps)
  names(overlaps) <- onames
  overlaps$rel_frame <- getRelativeFrames(overlaps)
  rel_frame <- getRelativeFrames(overlaps)
  
  overlaps <- subsetByOverlaps(overlaps, display_range)
  
  intersections <- trimOverlaps(overlaps,display_range)
  intersections <- groupGRangesBy(intersections)
  
  locations <- pmapToTranscriptF(intersections, display_range)
  layers <- geneTrackLayer(locations)
  
  locations <- unlistGrl(locations)
  ###
  locations <- IRanges::mid(locations)
  rel_locations <- ((locations / plot_width) * 0.95) - 0.475
  
  
  p <- ggplot(ct, aes(x = library, group = library, y  = log2(RPKM))) +
    geom_violin() +
    geom_boxplot(width=0.1) +
    geom_point(data = ct[transcript %in% names(overlaps)], mapping = aes(x = library, fill = transcript, color = transcript, y  = log2(RPKM),text = transcript), position = position_nudge(x = rel_locations), size = 4, shape = 15) +
    theme_bw() 
    # theme(legend.title = element_text(size = 16),
    # legend.text = element_text(size = 14) ) 

  p <- ggplotly(p)
  p <- p  %>% plotly::config(
    toImageButtonOptions = list(
      format = "svg",
      filename = "distribution_plot"))
  p %>% layout(xaxis = list(title = list(font = list(size = 40)), tickfont = list(size = 22)),
         yaxis = list(title = list(font = list(size = 40)), tickfont = list(size = 22)))
  return(ggplotly(p))
  
}
