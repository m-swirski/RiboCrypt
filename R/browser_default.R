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
#' @param kmers numeric (integer), bin positions into kmers. Default NULL, which is equal to 1,
#' i.e. no binning.
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
#' @param summary_track_type character, default is same as 'frames_type'
#'  argument
#' @param annotation_proportions numeric vector of length 3 (seq displayed), or 2 (seq not displayed). Relative sizes of annotation tracks.
#' @param summary_track logical, default FALSE. Display a top track, that is the sum
#' of all tracks.
#' @param AA_code Genetic code for amino acid display. Default is SGC0 (standard: Vertebrate).
#' See \code{Biostrings::GENETIC_CODE_TABLE} for options. To change to bacterial, do:
#' \code{Biostrings::getGeneticCode("11")}
#' @param log_scale logical, default FALSE. Log2 scale the count values, for easier visualization of shapes.
#' @param BPPARAM how many cores/threads to use? default: \code{BiocParallel::SerialParam()}.
#'  To see number of threads used for multicores, do \code{BiocParallel::bpparam()$workers}.
#'  You can also add a time remaining bar, for a more detailed pipeline.
#' @param export.format character, default: "svg". alternative: "png".
#' when you click the top right image button export, what should it export as?
#' @param frames_subset character, default "all". Alternatives: "red", "green", "blue".
#' @param zoom_range numeric or NULL. Default NULL, else an interval of xrange of the plot to
#'  initate plot on and highlight.
#' @param tx_annotation NULL
#' @param collapse_intron_flank integer, default 100
#' @param frame_colors character, color theme of the 3 coding frames, default "R", else "Color_blind".
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
#' # With Coding frame coloring
#' multiOmicsPlot_list(mrna[1], annotation = cds[1], reference_sequence = findFa(df),
#'                     frames_type = "columns", leader_extension = 30, trailer_extension = 30,
#'                     withFrames = TRUE,
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
                                log_scale = FALSE, BPPARAM = BiocParallel::SerialParam(), summary_track = FALSE,
                                summary_track_type = frames_type,
                                export.format = "svg", frames_subset = "all",
                                zoom_range = NULL, tx_annotation = NULL, collapse_intron_flank = 100,
                                frame_colors = "R") {

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
    custom_motif, log_scale, BPPARAM, "",
    summary_track, summary_track_type, export.format, frames_subset,
    zoom_range, tx_annotation, collapse_intron_flank, frame_colors)

}
