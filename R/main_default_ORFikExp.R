
#' Multi-omics plot using ORFik experiment input
#'
#' Customizable html plots for visualizing genomic data.
#' @inheritParams multiOmicsPlot_list
#' @param df an ORFik \code{\link[ORFik]{experiment}} or a list containing ORFik experiments.
#' Usually a list when you have split Ribo-seq and RNA-seq etc.
#' @param reference_sequence the genome reference, default ORFik::findFa(df)
#' @param reads the NGS libraries, as a list of \code{\link[GenomicRanges]{GRanges}}
#' with or without 'score' column for replicates. Can also be a covRle object
#' of precomputed coverage.
#' Default: \code{outputLibs(df, type = "pshifted", output.mode = "envirlist",
#' naming = "full", BPPARAM = BiocParallel::SerialParam())}
#' @param withFrames a logical vector, default
#' \code{libraryTypes(df, uniqueTypes = FALSE) \%in\% c("RFP", "RPF", "LSU")}
#' Alternative: a length 1 or same length as list length of "reads" argument.
#' @param ylabels character, default \code{bamVarName(df)}. Name of libraries in "reads" list argument.
#' @param plot_name character, default "default" (will create name from display_range name).
#' @param input_id character path, default: "", id for shiny to disply structures,
#'  should be "" for local users.
#' @return the plot object
#' @importFrom GenomicFeatures extractTranscriptSeqs seqlevels<-
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom GenomicRanges GRangesList mcols
#' @importFrom htmlwidgets onRender
#' @export
#' @examples
#' library(RiboCrypt)
#' df <- ORFik.template.experiment()[9,] #Use third library in experiment only
#' cds <- loadRegion(df, "cds")
#' multiOmicsPlot_ORFikExp(extendLeaders(extendTrailers(cds[1], 30), 30), df = df,
#'                         frames_type = "columns")
multiOmicsPlot_ORFikExp <- function(display_range, df, annotation = "cds",reference_sequence = findFa(df),
                                reads = outputLibs(df, type = "pshifted", output.mode = "envirlist",
                                                   naming = "full", BPPARAM = BiocParallel::SerialParam()),
                                viewMode = c("tx", "genomic")[1],
                                custom_regions = NULL,
                                leader_extension = 0, trailer_extension = 0,
                                withFrames = libraryTypes(df, uniqueTypes = FALSE) %in% c("RFP", "RPF", "LSU"),
                                frames_type = "lines", colors = NULL, kmers = NULL, kmers_type = c("mean", "sum")[1],
                                ylabels = bamVarName(df), lib_to_annotation_proportions = c(0.8,0.2),lib_proportions = NULL,
                                annotation_proportions = NULL,width = NULL, height = NULL,
                                plot_name = "default", plot_title = NULL,
                                display_sequence = c("both","nt", "aa", "none")[1], seq_render_dist = 100,
                                aa_letter_code = c("one_letter", "three_letters")[1],
                                annotation_names = NULL, start_codons = "ATG", stop_codons = c("TAA", "TAG", "TGA"),
                                custom_motif = NULL, BPPARAM = BiocParallel::SerialParam(),
                                input_id = "", summary_track = FALSE) {

  multiOmicsPlot_internal(display_range, df, annotation,reference_sequence,
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
                          custom_motif, BPPARAM, input_id, summary_track)

}
