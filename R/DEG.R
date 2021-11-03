#' Differential expression plots (1D or 2D)
#'
#' Gives you interactive 1D or 2D DE plots
#' @param dt a data.table with results from a differential
#' expression run. Normally from: \code{ORFik::DTEG.analysis(df1, df2)}
#' @param draw_non_regulated logical, default FALSE.
#' Should non-regulated rows be included in the plot?
#' Will make the plot faster to render if skipped (FALSE)
#' @return plotly object
#' @export
#' @examples
#' df <- ORFik.template.experiment()
#' dt <- DTEG.analysis(df[df$libtype == "RFP",], df[df$libtype == "RNA",],
#'                     output.dir = NULL)
#' #DEG_plot(dt, draw_non_regulated = TRUE)
DEG_plot <- function(dt, draw_non_regulated = FALSE) {
  stopifnot(is(dt, "data.table"))
  columns_must_exists <- c("variable", "Regulation", "id", "rna", "rfp", "te")
  if (!(all(columns_must_exists %in% colnames(dt))))
    stop("dt must minimally contain the columns: ",
         paste(columns_must_exists, collapse = ", "))
  if (!draw_non_regulated) dt <- dt[Regulation != "No change",]
  if (nrow(dt) == 0) {
    warning("dt input had no valid rows to plot")
    return(invisible(NULL))
  }
  id <- NULL
  dt <- highlight_key(dt, ~id, "Select a transcript")
  gg <- ggplot(dt, aes(x = rna, y = rfp, color = Regulation, frame = variable, id = id)) +
    geom_vline(aes(xintercept = 0), color = "red")+
    geom_hline(aes(yintercept = 0), color = "red") +
    geom_abline(slope = 1, intercept = 0, color = "grey", linetype = 2)+
    geom_point() +
    xlim(c(-1,1))+
    ylim(c(-1,1))

  suppressWarnings(select <- highlight(
    ggplotly(gg, tooltip = c("id")) %>% toWebGL() %>% partial_bundle(),
    on = "plotly_selected",
    selectize = TRUE, persistent = TRUE
  ))
  return(select)
}
