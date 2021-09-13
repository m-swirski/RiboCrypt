#' Multi-omics plot using ORFik experiment input
#'
#' Customizable html plots for visualizing genomic data.
#' @inheritParams multiOmicsPlot_list
#' @param df an ORFik \code{\link{experiment}} or a list containing ORFik experiments.
#' Usually a list when you have split Ribo-seq and RNA-seq etc.
#' @param reference_sequence the genome reference, default findFa(df)
#' @param reads the NGS libraries, as a list of \code{\link{GRanges}} with or without score column for replicates.
#' Default: \code{outputLibs(df, type = "pshifted", output.mode = "envirlist", naming = "full")}
#' @param withFrames a logical vector, default
#' \code{libraryTypes(df, uniqueTypes = FALSE) \%in\% c("RFP", "RPF", "LSU")}
#' Alternative: a length 1 or same length as list length of "reads" argument.
#' @param frames_type character, default "lines". Alternative:\cr
#' - columns \cr
#' - stacks \cr
#' - area
#' @param colors character, default NULL (automatic colouring). If "withFrames" argument
#' is TRUE, colors are set to to c("red", "green", "blue") for the 3 frames.
#' Alternative:  Character vector of length 1 or length of "reads" list argument.
#' @param kmers numeric (integer), bin positions into kmers.
#' @param ylabels character, default \code{bamVarName(df)}. Name of libraries in "reads" list arugment.
#' @param proportions numeric, default NULL. Width of plot.
#' @param width numeric, default NULL. Width of plot.
#' @param height numeric, default NULL. Height of plot.
#' @param plot_name = character, default "default" (will create name from target_range name).
#' Alternative: custom name for region.
#' @param plot_title character, default NULL. A title for plot.
#' @param display_sequence logical, default FALSE. If TRUE, display nucleotide sequence in plot.
#' @param annotation_names character, default NULL. Alternative naming for annotation.
#' @return the plot object
#' @importFrom GenomicFeatures extractTranscriptSeqs
#' @export
multiOmicsPlot_ORFikExp <- function(target_range, annotation = target_range, df, reference_sequence = findFa(df),
                                    reads = outputLibs(df, type = "pshifted", output.mode = "envirlist", naming = "full"),
                                    withFrames = libraryTypes(df, uniqueTypes = FALSE) %in% c("RFP", "RPF", "LSU"),
                                    frames_type = "lines", colors = NULL, kmers = NULL,
                                    ylabels = bamVarName(df), proportions = NULL, width = NULL, height = NULL,
                                    plot_name = "default", plot_title = NULL, display_sequence = FALSE,
                                    annotation_names = NULL) {
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
