#' Multi-omics plot using ORFik experiment input
#'
#' Customizable html plots for visualizing genomic data.
#' @inheritParams multiOmicsPlot_list
#' @param df an ORFik \code{\link[ORFik]{experiment}} or a list containing ORFik experiments.
#' Usually a list when you have split Ribo-seq and RNA-seq etc.
#' @param reference_sequence the genome reference, default ORFik::findFa(df)
#' @param reads the NGS libraries, as a list of \code{\link[GenomicRanges]{GRanges}} with or without score column for replicates.
#' Default: \code{outputLibs(df, type = "pshifted", output.mode = "envirlist", naming = "full")}
#' @param withFrames a logical vector, default
#' \code{libraryTypes(df, uniqueTypes = FALSE) \%in\% c("RFP", "RPF", "LSU")}
#' Alternative: a length 1 or same length as list length of "reads" argument.
#' @param ylabels character, default \code{bamVarName(df)}. Name of libraries in "reads" list argument.
#' @param plot_name character, default "default" (will create name from target_range name).
#' Alternative: custom name for region.
#' @return the plot object
#' @importFrom GenomicFeatures extractTranscriptSeqs seqlevels<-
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom GenomicRanges GRangesList mcols
#' @export
#' @examples
#' library(ORFik)
#' df <- ORFik.template.experiment()[3,] #Use third library in experiment only
#' if (requireNamespace("BSgenome.Hsapiens.UCSC.hg19")) {
#'   cds <- loadRegion(df, "cds")
#'   multiOmicsPlot_ORFikExp(extendLeaders(extendTrailers(cds[1], 30), 30), df = df,
#'                         reference_sequence = BSgenome.Hsapiens.UCSC.hg19::Hsapiens,
#'                         frames_type = "columns")
#' }
multiOmicsPlot_ORFikExp <- function(target_range, annotation = target_range, df, reference_sequence = findFa(df),
                                    reads = outputLibs(df, type = "pshifted", output.mode = "envirlist", naming = "full"),
                                    withFrames = libraryTypes(df, uniqueTypes = FALSE) %in% c("RFP", "RPF", "LSU"),
                                    frames_type = "lines", colors = NULL, kmers = NULL, kmers_type = c("mean", "sum")[1],
                                    ylabels = bamVarName(df), proportions = NULL, width = NULL, height = NULL,
                                    plot_name = "default", plot_title = NULL, display_sequence = FALSE,
                                    annotation_names = NULL, start_codons = "ATG", stop_codons = c("TAA", "TAG", "TGA"),
                                    custom_motif = NULL) {
  multiOmicsController()
  target_seq <- extractTranscriptSeqs(reference_sequence, target_range)
  seq_panel <- createSeqPanel(target_seq[[1]], start_codons = start_codons,
                              stop_codons = stop_codons, custom_motif = custom_motif)



  gene_model_panel <- createGeneModelPanel(target_range, annotation)
  lines <- gene_model_panel[[2]]
  gene_model_panel <- gene_model_panel[[1]]
  plots <- mapply(function(x,y,z,c,g) createSinglePlot(target_range, x,y,z,c,kmers_type, g, lines, type = frames_type), reads, withFrames, colors, kmers, ylabels, SIMPLIFY = FALSE)

  if (!display_sequence){
    plots <- c(plots, list(automateTicksGMP(gene_model_panel), automateTicksX(seq_panel)))
    multiomics_plot <- subplot(plots,
                               margin = 0,
                               nrows = length(reads) + 2,
                               heights = proportions,
                               shareX = TRUE,
                               titleY = TRUE,
                               titleX = TRUE)
  } else {
    letters <- nt_bar(target_seq)
    plots <- c(plots, list(automateTicksLetters(letters),automateTicksGMP(gene_model_panel), automateTicksX(seq_panel)))
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
