
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
#' @param plot_name character, default "default" (will create name from display_range name).
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
                                reads = outputLibs(df, type = "pshifted", output.mode = "envirlist", naming = "full"),
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
                                custom_motif = NULL, BPPARAM = BiocParallel::SerialParam()) {
  multiOmicsController()
  # Get sequence and create basic seq panel
  target_seq <- extractTranscriptSeqs(reference_sequence, display_range)
  seq_panel <- createSeqPanel(target_seq[[1]], start_codons = start_codons,
                              stop_codons = stop_codons, custom_motif = custom_motif)


  # Get the panel for the annotation track
  gene_model_panel <- createGeneModelPanel(display_range, annotation,
                                           custom_regions = custom_regions,
                                           viewMode = viewMode)
  lines <- gene_model_panel[[2]]
  gene_model_panel <- gene_model_panel[[1]]
  # Get NGS data tracks
  force(display_range)
  force(kmers_type)
  force(lines)
  force(frames_type)
  force(reads)
  force(withFrames)
  force(colors)
  force(kmers)
  force(ylabels)
  if (is(BPPARAM, "SerialParam")) {
    plots <- mapply(function(x,y,z,c,g) createSinglePlot(display_range, x,y,z,c,kmers_type, g, lines, type = frames_type),
                      reads, withFrames, colors, kmers, ylabels, SIMPLIFY = FALSE)
  } else {
    plots <- bpmapply(function(x,y,z,c,g) createSinglePlot(display_range, x,y,z,c,kmers_type, g, lines, type = frames_type),
                      reads, withFrames, colors, kmers, ylabels, SIMPLIFY = FALSE, BPPARAM = BPPARAM)
  }


  if (display_sequence %in% c("none", FALSE)) { # plotly subplot without sequence track
    plots <- c(plots, list(automateTicksGMP(gene_model_panel), automateTicksX(seq_panel)))
    multiomics_plot <- subplot(plots,
                               margin = 0,
                               nrows = length(reads) + 2,
                               heights = proportions,
                               shareX = TRUE,
                               titleY = TRUE,
                               titleX = TRUE)
  } else { # plotly subplot with sequence track
    nplots <- length(plots)
    nt_area <- ggplot() +
      theme(axis.title = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank()) +
      theme(plot.margin = unit(c(0,0,0,0), "pt"))+
      scale_x_continuous(expand = c(0,0))

    plots <- c(plots, list(automateTicks(nt_area), automateTicksGMP(gene_model_panel),
                           automateTicksX(seq_panel)))
    multiomics_plot <- subplot(plots,
                               margin = 0,
                               nrows = length(reads) + 3,
                               heights = proportions,
                               shareX = TRUE,
                               titleY = TRUE,
                               titleX = TRUE)
    # Create sequence zoom logic (javascript)
    display_dist <- nchar(target_seq)
    render_on_zoom_data <- fetch_JS_seq(target_seq = target_seq, nplots = nplots,
                            distance = seq_render_dist, display_dist = display_dist, aa_letter_code = aa_letter_code)
    select_region_on_click_data <- list(nplots = nplots)
    multiomics_plot <- multiomics_plot %>%
      onRender(RiboCrypt:::fetchJS("render_on_zoom.js"), render_on_zoom_data) %>%
      onRender(RiboCrypt:::fetchJS("select_region_on_click.js"), select_region_on_click_data)
  }

  multiomics_plot <- multiomics_plot %>% plotly::config(
    toImageButtonOptions = list(
      format = "svg",
      filename = ifelse(plot_name == "default", names(display_range), plot_name),
      width = width,
      height = height))
  if (!is.null(plot_title)) multiomics_plot <- multiomics_plot %>% plotly::layout(title = plot_title)

  return(multiomics_plot)
}
