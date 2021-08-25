library(GenomicFeatures)
library(plotly)
library(ORFik)
library(ggplot2)
library(data.table)
library(GenomicAlignments)
library(BSgenome)

multiOmicsPlot_list <- function(target_range, annotation, reference_sequence, reads, withFrames = NULL, frames_type = "lines", colors = NULL, kmers = NULL, ylabels = NULL, proportions = NULL, width = NULL, height = NULL,plot_name = "default", plot_title = NULL, display_sequence = FALSE, annotation_names = NULL) {
  seqlevels(target_range) <- seqlevels(annotation)
  target_range <- GRangesList(target_range)
  
  # if (is(annotation, "GRangesList")) annotation <- unlist(annotation)
  
  if (!is.null(annotation_names)) {
    if (length(annotation_names) == 1){
      if (annotation_names %in% annotation) {
        names(annotation) <- mcols(annotation)[[annotation_names]]
      } else stop(wmsg("wrong annotation_names argument"))
    } else if (length(annotation_names) == length(annotation)) {
      names(annotation) <- annotation_names
    } else stop(wmsg("wrong annotation_names argument"))
  }
  
  if(class(reads) != "list") reads <- list(reads)
  
  if (length(withFrames) == 0) withFrames <- FALSE
  if (!(length(withFrames)  %in% c(1, length(reads)))) stop("length of withFrames must be 0, 1 or the same as reads list")
  if (length(withFrames) == 1) withFrames <- rep(withFrames, length(reads))
  
  if (length(colors) == 0) colors <- 1:length(reads)
  if (!(length(colors)  %in% c(1, length(reads)))) stop("length of colors must be 0, 1 or the same as reads list")
  if (length(colors) == 1) colors <- rep(colors, length(reads))
  
  if (length(kmers) == 0) kmers <- 1
  if (!(length(kmers)  %in% c(1, length(reads)))) stop("length of kmers must be 0, 1 or the same as reads list")
  if (length(kmers) == 1) kmers <- rep(kmers, length(reads))
  
  if (length(ylabels) == 0) ylabels <- as.character(1:length(reads))
  if (!(length(ylabels)  %in% c(1, length(reads)))) stop("length of ylabels must be 0, 1 or the same as reads list")
  if (length(ylabels) == 1) ylabels <- rep(ylabels, length(reads))
  
  if (length(proportions) == 0) proportions <- 1
  if (!(length(proportions)  %in% c(1, length(reads)))) stop("length of proportions must be 0, 1 or the same as reads list")
  if (length(proportions) == 1) proportions <- rep(proportions, length(reads))
  
  if (!display_sequence) {
    max_sum <- 1  - sum(0.03,0.1)
    proportions <- proportions / (sum(proportions) / max_sum)
    proportions <- c(proportions, c(0.03,0.1))
  } else {
    max_sum <- 1  - sum(0.035, 0.03,0.1)
    proportions <- proportions / (sum(proportions) / max_sum)
    proportions <- c(proportions, c(0.035, 0.03,0.1))
  }
  
  
  target_seq <- extractTranscriptSeqs(reference_sequence, target_range)
  seq_panel <- createSeqPanel(target_seq[[1]])
  
  
  
  gene_model_panel <- createGeneModelPanel(target_range, annotation)
  lines <- gene_model_panel[[2]]
  gene_model_panel <- gene_model_panel[[1]]
  plots <- mapply(function(x,y,z,c,g) createSinglePlot(target_range, x,y,z,c,g, lines, type = frames_type), reads, withFrames, colors, kmers, ylabels, SIMPLIFY = FALSE)
  
  if (!display_sequence){ 
    plots <- c(plots, list(automateTicks(gene_model_panel), automateTicksX(seq_panel)))
    multiomics_plot <- subplot(plots, 
                               margin = 0,
                               nrows = length(reads) + 2, 
                               heights = proportions, 
                               shareX = TRUE,
                               titleY = TRUE, 
                               titleX = TRUE)
  } else {
    letters <- nt_bar(target_seq)
    plots <- c(plots, list(automateTicks(letters),automateTicks(gene_model_panel), automateTicksX(seq_panel)))
    multiomics_plot <- subplot(plots,
                               margin = 0,
                               nrows = length(reads) + 3, 
                               heights = proportions, 
                               shareX = TRUE,
                               titleY = TRUE, 
                               titleX = TRUE)
  }
  
  multiomics_plot <- multiomics_plot %>% plotly::config(
    toImageButtonOptions = list(
      format = "svg",
      filename = ifelse(plot_name == "default",names(target_range),plot_name),
      width = width,
      height = height))
  if (!is.null(plot_title)) multiomics_plot <- multiomics_plot %>% plotly::layout(title = plot_title)
  return(multiomics_plot)
  
}


multiOmicsPlot_animate <- function(target_range, annotation, reference_sequence, reads, withFrames = NULL, colors = NULL, kmers = NULL, ylabels = NULL, proportions = NULL, width = NULL, height = NULL,plot_name = "default", plot_title = NULL, display_sequence = FALSE, annotation_names = NULL) {
  seqlevels(target_range) <- seqlevels(annotation)
  target_range <- GRangesList(target_range)
  
  
  # if (is(annotation, "GRangesList")) annotation <- unlist(annotation)
  if (!is.null(annotation_names)) {
    if (length(annotation_names) == 1){
      if (annotation_names %in% annotation) {
        names(annotation) <- mcols(annotation)[[annotation_names]]
      } else stop(wmsg("wrong annotation_names argument"))
    } else if (length(annotation_names) == length(annotation)) {
      names(annotation) <- annotation_names
    } else stop(wmsg("wrong annotation_names argument"))
  }
  
  if(class(reads) != "list") reads <- list(reads)
  
  if (length(withFrames) == 0) withFrames <- FALSE
  if (!(length(withFrames)  %in% c(1, length(reads)))) stop("length of withFrames must be 0, 1 or the same as reads list")
  if (length(withFrames) == 1) withFrames <- rep(withFrames, length(reads))
  
  if (length(colors) == 0) colors <- 1:length(reads)
  if (!(length(colors)  %in% c(1, length(reads)))) stop("length of colors must be 0, 1 or the same as reads list")
  if (length(colors) == 1) colors <- rep(colors, length(reads))
  
  if (length(kmers) == 0) kmers <- 1
  if (!(length(kmers)  %in% c(1, length(reads)))) stop("length of kmers must be 0, 1 or the same as reads list")
  if (length(kmers) == 1) kmers <- rep(kmers, length(reads))
  
  if (length(ylabels) == 0) ylabels <- as.character(1:length(reads))
  if (!(length(ylabels)  %in% c(1, length(reads)))) stop("length of ylabels must be 0, 1 or the same as reads list")
  if (length(ylabels) == 1) ylabels <- rep(ylabels, length(reads))
  
  if (length(proportions) == 0) proportions <- 1
  if (!(length(proportions)  %in% c(1, length(reads)))) stop("length of proportions must be 0, 1 or the same as reads list")
  if (length(proportions) == 1) proportions <- rep(proportions, length(reads))
  
  if (!display_sequence) {
    max_sum <- 1  - sum(0.03,0.1)
    proportions <- proportions / (sum(proportions) / max_sum)
    proportions <- c(proportions, c(0.03,0.1))
  } else {
    max_sum <- 1  - sum(0.035, 0.03,0.1)
    proportions <- proportions / (sum(proportions) / max_sum)
    proportions <- c(proportions, c(0.035, 0.03,0.1))
  }
  
  
  target_seq <- extractTranscriptSeqs(reference_sequence, target_range)
  seq_panel <- createSeqPanel(target_seq[[1]], frame=1:length(reads))
  
  
  gene_model_panel <- createGeneModelPanel(target_range, annotation, frame = 1:length(reads))
  lines <- gene_model_panel[[2]]
  gene_model_panel <- gene_model_panel[[1]]
  
  profiles <- mapply(function(x,y,z) getProfileAnimate(target_range, x, y, z), reads, withFrames, kmers,  SIMPLIFY = FALSE)

  profiles <- rbindlist(profiles, idcol = "file")

  plot <- getPlotAnimate(profiles, withFrames = withFrames[1], colors = colors[1], ylabels = ylabels[1], lines = lines)

  plots <- list(plot, automateTicks(gene_model_panel), automateTicksX(seq_panel))
  
  if (!display_sequence){ 
    multiomics_plot <- subplot(plots, 
                               margin = 0,
                               nrows = 3, 
                               heights = c(0.8,0.05, 0.15), 
                               shareX = TRUE,
                               titleY = TRUE, 
                               titleX = TRUE)
    multiomics_plot <- lineDeSimplify(multiomics_plot)
  } else {
    letters <- nt_bar(target_seq)
    multiomics_plot <- subplot(c(plots,
                                 automateTicks(letters)),
                               margin = 0,
                               nrows = length(reads) + 3, 
                               heights = proportions, 
                               shareX = TRUE,
                               titleY = TRUE, 
                               titleX = TRUE)
    multiomics_plot <- lineDeSimplify(multiomics_plot)
    
  }
  
  multiomics_plot <- multiomics_plot %>% plotly::config(
    toImageButtonOptions = list(
      format = "svg",
      filename = ifelse(plot_name == "default",names(target_range),plot_name),
      width = width,
      height = height))
  if (!is.null(plot_title)) multiomics_plot <- multiomics_plot %>% plotly::layout(title = plot_title)
  return(multiomics_plot)
  
}
