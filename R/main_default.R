#' Multi-omics plot using list input
#'
#' Customizable html plots for visualizing genomic data.
#' @param display_range the whole region to visualize,
#'  a \code{\link[GenomicRanges]{GRangesList}} or \code{\link[GenomicRanges]{GRanges}} object
#' @param annotation the whole annotation which your target region is a subset,
#' a \code{\link[GenomicRanges]{GRangesList}} or \code{\link[GenomicRanges]{GRanges}} object
#' @param reference_sequence the genome reference,
#' a \code{\link[Rsamtools]{FaFile}} or \code{\link[Rsamtools]{FaFile}} convertible object
#' @param reads the NGS libraries, as a list of \code{\link[GenomicRanges]{GRanges}} with or without score column for replicates.
#' @param viewMode character, default "tx" (transcript coordinates, first position is 1,
#' exons are merged into a single sequence)\cr
#' Alternative: "genomic" (genomic coordinates, first position is first position in
#' \code{display_range} argument. Introns are displayed).
#' @param custom_regions a GRangesList or NULL, default: NULL.
#'  The alternative annotation, like self defined uORFs etc. The vertical annotation bars will have
#'  a different color.
#' @param leader_extension integer, default 0. (How much to extend view upstream)
#' @param trailer_extension integer, default 0. (How much to extend view downstream)
#' @param withFrames a logical vector, default NULL. Alternative: a length 1 or same length as list length of "reads" argument.
#' @param frames_type character, default "lines". Alternative:\cr
#' - columns \cr
#' - stacks \cr
#' - area
#' @param colors character, default NULL (automatic colouring). If "withFrames" argument
#' is TRUE, colors are set to to c("red", "green", "blue") for the 3 frames.
#' Alternative:  Character vector of length 1 or length of "reads" list argument.
#' @param kmers numeric (integer), bin positions into kmers.
#' @param kmers_type character, function used for kmers sliding window. default: "mean", alternative: "sum"
#' @param ylabels character, default NULL. Name of libraries in "reads" list arugment.
#' @param width numeric, default NULL. Width of plot.
#' @param height numeric, default NULL. Height of plot.
#' @param plot_name = character, default "default" (will create name from display_range name).
#' Alternative: custom name for region.
#' @param plot_title character, default NULL. A title for plot.
#' @param display_sequence character/logical, default \code{c("both","nt", "aa", "none")[1]}.
#' If TRUE or "both", display nucleotide and aa sequence in plot.
#' @param annotation_names character, default NULL. Alternative naming for annotation.
#' @param seq_render_dist integer, default  100. The sequences will appear after zooming below this threshold.
#' @param aa_letter_code character, when set to "three_letters", three letter amino acid code is used. One letter by default.
#' @param lib_to_annotation_proportions numeric vector of length 2. relative sizes of profiles and annotation.
#' @param lib_proportions numeric vector of length equal to displayed libs. Relative sizes of profiles displayed
#' @param annotation_proportions numeric vector of length 3 (seq displayed), or 2 (seq not displayed). Relative sizes of annotation tracks.
#' @param AA_code Genetic code for amino acid display. Default is SGC0 (standard: Vertebrate).
#' See \code{Biostrings::GENETIC_CODE_TABLE} for options. To change to bacterial, do:
#' \code{Biostrings::getGeneticCode("11")}
#' @param BPPARAM how many cores/threads to use? default: \code{BiocParallel::SerialParam()}.
#'  To see number of threads used for multicores, do \code{BiocParallel::bpparam()$workers}.
#'  You can also add a time remaining bar, for a more detailed pipeline.
#' @inheritParams createSeqPanelPattern
#' @return the plot object
#' @importFrom GenomicFeatures extractTranscriptSeqs
#' @importFrom BiocParallel bpparam bpmapply
#' @importFrom Biostrings GENETIC_CODE
#' @export
#' @examples
#' library(RiboCrypt)
#' df <- ORFik.template.experiment()[9:10,]
#' cds <- loadRegion(df, "cds")
#' mrna <- loadRegion(df, "mrna")
#' multiOmicsPlot_list(mrna[1], annotation = cds[1], reference_sequence = findFa(df),
#'                     frames_type = "columns", leader_extension = 30, trailer_extension = 30,
#'                     reads = outputLibs(df, type = "pshifted", output.mode = "envirlist",
#'                                   naming = "full", BPPARAM = BiocParallel::SerialParam()))
multiOmicsPlot_list <- function(display_range, annotation = display_range, reference_sequence,
                                reads, viewMode = c("tx", "genomic")[1], custom_regions = NULL,
                                leader_extension = 0, trailer_extension = 0,
                                withFrames = NULL,
                                frames_type = "lines", colors = NULL,
                                kmers = NULL, kmers_type = c("mean", "sum")[1],
                                ylabels = NULL, lib_to_annotation_proportions = c(0.8,0.2),
                                lib_proportions = NULL, annotation_proportions = NULL,
                                width = NULL, height = NULL, plot_name = "default",
                                plot_title = NULL,
                                display_sequence = c("both","nt", "aa", "none")[1], seq_render_dist = 100,
                                aa_letter_code = c("one_letter", "three_letters")[1],
                                annotation_names = NULL,
                                start_codons = "ATG", stop_codons = c("TAA", "TAG", "TGA"),
                                custom_motif = NULL, AA_code = Biostrings::GENETIC_CODE,
                                BPPARAM = BiocParallel::SerialParam()) {

  multiOmicsPlot_internal(display_range, df = NULL, annotation,reference_sequence,
    reads,
    viewMode,
    custom_regions,
    leader_extension, trailer_extension,
    withFrames,
    frames_type, colors, kmers, kmers_type,
    ylabels, lib_to_annotation_proportions,lib_proportions,
    annotation_proportions,width, height,
    plot_name, plot_title,
    display_sequence, seq_render_dist,
    aa_letter_code,
    annotation_names, start_codons, stop_codons,
    custom_motif, BPPARAM)

}

#' Multi-omics animation using list input
#'
#' The animation will move with a play butten, there is
#' 1 transition per library given.
#' @inheritParams multiOmicsPlot_list
#' @return the plot object
#' @export
#' @examples
#' library(RiboCrypt)
#' df <- ORFik.template.experiment()[9:10,]
#' cds <- loadRegion(df, "cds")
#' mrna <- loadRegion(df, "mrna")
#' # multiOmicsPlot_animate(mrna[1], annotation = cds[1], reference_sequence = findFa(df),
#' #                     frames_type = "columns", leader_extension = 30, trailer_extension = 30,
#' #                     reads = outputLibs(df, type = "pshifted", output.mode = "envirlist",
#' #                                   naming = "full", BPPARAM = BiocParallel::SerialParam()))
multiOmicsPlot_animate <- function(display_range, annotation = display_range, reference_sequence,
                                   reads, viewMode = c("tx", "genomic")[1], custom_regions = NULL,
                                   leader_extension = 0, trailer_extension = 0,
                                   withFrames = NULL,
                                   frames_type = "lines", colors = NULL,
                                   kmers = NULL, kmers_type = c("mean", "sum")[1],
                                   ylabels = NULL, lib_to_annotation_proportions = c(0.8,0.2),
                                   lib_proportions = NULL, annotation_proportions = NULL,
                                   width = NULL, height = NULL, plot_name = "default",
                                   plot_title = NULL,
                                   display_sequence = c("both","nt", "aa", "none")[1], seq_render_dist = 100,
                                   aa_letter_code = c("one_letter", "three_letters")[1],
                                   annotation_names = NULL,
                                   start_codons = "ATG", stop_codons = c("TAA", "TAG", "TGA"),
                                   custom_motif = NULL,
                                   BPPARAM = BiocParallel::SerialParam()) {
  multiOmicsController()
  # Get sequence and create basic seq panel
  target_seq <- extractTranscriptSeqs(reference_sequence, display_range)
  read_names <- names(reads) # Names for frames for seq panels
  seq_panel <- createSeqPanelPattern(
    target_seq[[1]], start_codons = start_codons,
    stop_codons = stop_codons, custom_motif = custom_motif,
    frame = read_names)


  # Get the panel for the annotation track
  gene_model_panel <- createGeneModelPanel(display_range, annotation,
                                           custom_regions = custom_regions,
                                           viewMode = viewMode,
                                           frame = read_names)
  lines <- gene_model_panel[[2]]
  gene_model_panel <- gene_model_panel[[1]]

  # profiles <- mapply(function(x,y,z) getProfileAnimate(display_range, x, y, z, kmers_type),
  #                    reads, withFrames, kmers,  SIMPLIFY = FALSE)
  profiles <- bpmapply(function(x,y,z) getProfileAnimate(display_range, x, y, z, kmers_type),
                        reads, withFrames, kmers,  SIMPLIFY = FALSE, BPPARAM = BPPARAM)

  profiles <- rbindlist(profiles, idcol = "file")

  plot <- list(getPlotAnimate(profiles, withFrames = withFrames[1],
                         colors = colors[1], ylabels = ylabels[1], lines = lines))

  if (display_sequence %in% c("none", FALSE)){
    plots <- c(plot, list(automateTicksGMP(gene_model_panel), automateTicksX(seq_panel)))
    multiomics_plot <- subplot(plots,
                               margin = 0,
                               nrows = 3,
                               heights = c(0.8, 0.05, 0.15),
                               shareX = TRUE,
                               titleY = TRUE, titleX = TRUE)
  } else {
    nplots <- length(plot)
    nt_area <- ggplot() +
      theme(axis.title = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank()) +
      theme(plot.margin = unit(c(0,0,0,0), "pt"))+
      scale_x_continuous(expand = c(0,0))
    plots <- c(plot, list(automateTicks(nt_area), automateTicksGMP(gene_model_panel),
                           automateTicksX(seq_panel)))
    multiomics_plot <- subplot(plots,
                               margin = 0,
                               nrows = 4,
                               heights = c(0.7, 0.05,  0.15, 0.1),
                               shareX = TRUE,
                               titleY = TRUE,
                               titleX = TRUE)
  }
  multiomics_plot <- lineDeSimplify(multiomics_plot)
  multiomics_plot <- multiomics_plot %>% plotly::config(
    toImageButtonOptions = list(
      format = "svg",
      filename = ifelse(plot_name == "default",names(display_range),plot_name),
      width = width,
      height = height))
  if (!is.null(plot_title)) multiomics_plot <- multiomics_plot %>% plotly::layout(title = plot_title)
  return(multiomics_plot)
}
