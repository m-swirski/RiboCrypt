#' Differential expression plots (1D or 2D)
#'
#' Gives you interactive 1D or 2D DE plots
#' @param dt a data.table with results from a differential
#' expression run. Normally from: \code{ORFik::DTEG.analysis(df1, df2)}
#' @param draw_non_regulated logical, default FALSE.
#' Should non-regulated rows be included in the plot?
#' Will make the plot faster to render if skipped (FALSE)
#' @param two_dimensions logical, default:
#' \code{ifelse("LFC" \%in\% colnames(dt), FALSE, TRUE)}
#' Is this two dimensional, like: Ribo-seq vs RNA-seq. Alternative, FALSE: Then
#' Log fold change vs mean counts
#' @param xlim numeric vector or character preset, default:
#'  \code{ifelse(two_dimensions, "bidir.max", "auto")}
#'  (Equal in both + / - direction, using max value + 0.5 of
#'  meanCounts(in 1d) / rna(in 2d) column of dt).
#'  If you want ggplot to decide limit, set to "auto". For numeric vector,
#'  specify min and max x limit: like c(-5, 5)
#' @param ylim numeric vector or character preset, default: "bidir.max"
#'  (Equal in both + / - direction, using max value + 0.5 of
#'  LFC(in 1d) / rfp(in 2d) column of dt).
#'  If you want ggplot to decide limit, set to "auto". For numeric vector,
#'  specify min and max x limit: like c(-5, 5)
#' @param xlab character, default:
#' \code{ifelse(two_dimensions, "RNA fold change (log2)", "Mean counts (log2)")}
#' @param ylab character, default:
#' \code{ifelse(two_dimensions, "RFP fold change (log2)",  "Fold change (log2)")}
#' @param color.values named character vector, default: \code{c("No change" = "black", "Significant" = "red",
#'  "Buffering" = "purple", "mRNA abundance" = "darkgreen",
#'  "Expression" = "blue", "Forwarded" = "yellow",
#'  "Inverse" = "aquamarine", "Translation" = "orange4")}
#' @return plotly object
#' @export
#' @examples
#' # Load experiment
#' df <- ORFik.template.experiment()
#' # 1 Dimensional analysis
#' dt <- DEG.analysis(df[df$libtype == "RNA",],
#'                     output.dir = NULL)
#' dt$Regulation[1] <- "Significant" # Fake sig level
#' DEG_plot(dt, draw_non_regulated = TRUE)
#' # 2 Dimensional analysis
#' dt_2d <- DTEG.analysis(df[df$libtype == "RFP",], df[df$libtype == "RNA",],
#'                     output.dir = NULL)
#' dt_2d$Regulation[4] <- "Translation" # Fake sig level
#' dt_2d$Regulation[5] <- "Buffering" # Fake sig level
#' DEG_plot(dt_2d, draw_non_regulated = TRUE)
DEG_plot <- function(dt, draw_non_regulated = FALSE,
                     xlim = ifelse(two_dimensions, "bidir.max", "auto"),
                     ylim = "bidir.max",
                     xlab = ifelse(two_dimensions, "RNA fold change (log2)", "Mean counts (log2)"),
                     ylab = ifelse(two_dimensions, "RFP fold change (log2)",  "Fold change (log2)"),
                     two_dimensions = ifelse("LFC" %in% colnames(dt), FALSE, TRUE),
                     color.values = c("No change" = "black", "Significant" = "red",
                                      "Buffering" = "purple", "mRNA abundance" = "darkgreen",
                                      "Expression" = "blue", "Forwarded" = "yellow",
                                      "Inverse" = "aquamarine", "Translation" = "orange4")) {
  stopifnot(is(dt, "data.table"))
  if (is.character(xlim)) stopifnot(xlim %in% c("bidir.max", "auto"))
  if (is.character(ylim)) stopifnot(ylim %in% c("bidir.max", "auto"))
  # Remove this line in next bioc
  if("variable" %in% colnames(dt)) colnames(dt) <- gsub("variable", "contrast", colnames(dt))
  columns_must_exists <- if (two_dimensions) {
    c("contrast", "Regulation", "id", "rna", "rfp", "te")
  } else c("contrast", "Regulation", "id", "LFC", "meanCounts")

  if (!(all(columns_must_exists %in% colnames(dt))))
    stop("dt must minimally contain the columns: ",
         paste(columns_must_exists, collapse = ", "))
  if (!draw_non_regulated) dt <- dt[Regulation != "No change",]
  if (nrow(dt) == 0) {
    warning("dt input had no valid rows to plot")
    return(invisible(NULL))
  }

  xlim <-
    if (!all(xlim == "auto")) {
      if (all(xlim == "bidir.max")) {
          if (two_dimensions) {c(-max(abs(dt$rna)) - 0.5, max(abs(dt$rna)) + 0.5)
          } else c(-max(abs(log2(dt$meanCounts))) - 0.5, max(log2(abs(dt$meanCounts))) + 0.5)
      } else xlim
    }
  ylim <-
    if (!all(ylim == "auto")) {
      if (all(ylim == "bidir.max")) {
          if (two_dimensions) {c(-max(abs(dt$rfp)) - 0.5, max(abs(dt$rfp)) + 0.5)
          } else c(-max(abs(dt$LFC)) - 0.5, max(abs(dt$LFC)) + 0.5)
      } else ylim
    }
  color.values <- color.values[as.character(unique(dt$Regulation))]
  id <- NULL
  dt <- highlight_key(dt, ~id, "Select a transcript")

  if (two_dimensions) {
    gg <- ggplot(dt, aes(x = rna, y = rfp, color = Regulation, frame = contrast, id = id)) +
      geom_vline(aes(xintercept = 0), color = "red") +
      geom_abline(slope = 1, intercept = 0, color = "grey", linetype = 2)
  } else { # One dimensional
    gg <- ggplot(dt, aes(x = log2(meanCounts), y = LFC, color = Regulation, frame = contrast, id = id))
  }
   gg <- gg +
    geom_hline(aes(yintercept = 0), color = "red") +
    geom_point() + xlab(xlab) + ylab(ylab) +
    scale_color_manual(values = color.values)
   if (!is.null(xlim)) gg <- gg + xlim(xlim)
   if (!is.null(ylim)) gg <- gg + ylim(ylim)

   hovertip <- "id"

   #vis <- ifelse(unique(dt$Regulation) == "", TRUE, "legendonly") #add_trace(visible = vis)
   select <- highlight(
     ggplotly(gg, tooltip = hovertip) %>%
       layout(autosize = TRUE) %>% toWebGL() %>% partial_bundle(),
     on = c('plotly_selected'), off = c('plotly_deselect'),
     selectize = TRUE, persistent = FALSE
   )
   return(select)
}
