DE_model_from_ctrl <- function(controls) {
  DE_model(controls()$dff, controls()$diff_method, controls()$full,
           controls()$target.contrast, controls()$design, controls()$all_libs,
           controls()$group_1, controls()$group_2)
}

DE_model <- function(df, method, other_tx, target.contrast = design[1],
                     design = ORFik::design(df), all_libs, group_1, group_2) {
  print("DE model activated!")
  if (nrow(df) < 2) stop("Differential expression only allowed for studies with > 1 sample")
  # TODO: add in design correctly
  if (any(group_1 %in% group_2))
    stop("Sample can not be in both contrast groups: ",
         paste(group_1[group_1 %in% group_2], collapse = " ,"))
  if (length(group_1) == 0) stop("Contrast group 1 is empty!")
  if (length(group_2) == 0) stop("Contrast group 2 is empty!")
  df <- df[all_libs %in% c(group_1, group_2),]
  if (nrow(df) < 2) stop("You selected < 2 samples, even though study has more!")
  counts <- countTable(df, "mrna", type = "summarized")
  counts <- counts[colnames(counts) %in% c(group_1, group_2)]
  if (!other_tx) {
    counts <- counts[filterTranscripts(df, 0, 1, 0)]
  }

  if (method == "DESeq2") {
    ORFik::DEG_model(df, counts = counts)
  } else if (method == "FPKM ratio") {
    ORFik::DEG_model_simple(df, target.contrast, design, counts = counts)
  } else stop("Invalid Differential method selected")
}

DE_model_results <- function(dt, controls, symbols_dt = symbols(controls()$dff)) {
  method <- controls()$diff_method
  if (method == "DESeq2") {
    dt <- DEG_model_results(dt,
                     controls()$target.contrast, pairs = controls()$pairs,
                     controls()$pval)
  } else if (method == "FPKM ratio") {
    print("Simplified DE model results")
    # TODO: Simplified model, make better
    regulation <- rep("No change", nrow(dt))
    regulation[abs(dt$LFC) * (dt$meanCounts + 1) > 100 & dt$meanCounts > 10 & abs(dt$LFC) > 0.5] <- "Significant"
    dt[, Regulation := regulation]
    dt[, LFC := round(LFC, 2)]
    dt[, meanCounts := round(meanCounts, 1)]
  } else stop("Invalid Differential method selected")
  return(append_gene_symbols(dt, symbols_dt))
}

#' Differential expression plots (1D or 2D)
#'
#' Gives you interactive 1D or 2D DE plots
#' @param dt a data.table with results from a differential
#' expression run. Normally from: \code{ORFik::DTEG.analysis(df1, df2)}
#' @param draw_non_regulated logical, default TRUE
#' Should non-regulated rows be included in the plot?
#' Will make the plot faster to render if skipped (FALSE)
#' @param add_search_bar logical, default TRUE. Add a crosstalk search bar
#' to search for genes in the plot
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
#' @param color.values named character vector, default: \code{c(
#'  "No change" = "black", "Significant" = "red",
#'  "Buffering" = "purple", "mRNA abundance" = "darkgreen",
#'  "Expression" = "blue", "Forwarded" = "yellow",
#'  "Inverse" = "aquamarine", "Translation" = "orange4")}
#' @param format character, default "png". Format for plotly bar.
#' @return plotly object or crosstalk bscols if add_search_bar is TRUE.
#' @importFrom crosstalk bscols SharedData filter_select
#' @export
#' @examples
#' # Load experiment
#' df <- ORFik.template.experiment()
#' df_rna <- df[df$libtype == "RNA",]
#' # 1 Dimensional analysis
#' dt <- DEG.analysis(df_rna)
#' dt$Regulation[1] <- "Significant" # Fake sig level
#' DEG_plot(dt, draw_non_regulated = TRUE)
#' # 2 Dimensional analysis
#' df_rfp <- df[df$libtype == "RFP",]
#' dt_2d <- DTEG.analysis(df_rfp, df_rna, output.dir = NULL)
#' dt_2d$Regulation[4] <- "Translation" # Fake sig level
#' dt_2d$rfp.lfc[4] <- -0.3 # Fake sig level
#' dt_2d$Regulation[5] <- "Buffering" # Fake sig level
#' dt_2d$rna.lfc[5] <- -0.3 # Fake sig level
#' DEG_plot(dt_2d, draw_non_regulated = TRUE)
#' # Add Gene symbols in ids for easier analysis
#' dt_2d_with_gene_ids <- ORFik::append_gene_symbols(dt_2d, symbols(df))
#' DEG_plot(dt_2d_with_gene_ids, draw_non_regulated = TRUE)
DEG_plot <- function(dt, draw_non_regulated = TRUE,
                     add_search_bar = TRUE,
                     xlim = ifelse(two_dimensions, "bidir.max", "auto"),
                     ylim = "bidir.max",
                     xlab = ifelse(two_dimensions, "RNA fold change (log2)", "Mean counts (log2)"),
                     ylab = ifelse(two_dimensions, "RFP fold change (log2)",  "Fold change (log2)"),
                     two_dimensions = ifelse("LFC" %in% colnames(dt), FALSE, TRUE),
                     color.values = c("No change" = "black", "Significant" = "red",
                                      "Buffering" = "purple", "mRNA abundance" = "darkgreen",
                                      "Expression" = "blue", "Forwarded" = "yellow",
                                      "Inverse" = "aquamarine", "Translation" = "orange4"),
                     format = "png") {
  DEG_plot_input_validation()
  if (nrow(dt) == 0) return(invisible(NULL))

  dt[, y_axis := if("rfp.lfc" %in% colnames(dt)) {rfp.lfc} else LFC]
  dt[, x_axis := if("rna.lfc" %in% colnames(dt)) {rna.lfc} else log2(meanCounts)]
  dt[, trace := paste0("id: ", id, "\nReg: ", Regulation)]


  if (two_dimensions) {
    values <- abs(c(dt$y_axis, dt$x_axis))
    max <- max(values[is.finite(values)], na.rm = TRUE) + 0.5
    max <- c(-max, max)
  } else {
    values <- dt$x_axis
    values <- values[is.finite(values)]
    max <- c(min(c(-0.5, values), na.rm = TRUE), max(values + 0.5, na.rm = TRUE))
  }
  return(DEG_plotly(dt, color.values, format, two_dimensions, max, xlab, ylab,
                    add_search_bar))
}

DEG_plotly <- function(dt, color.values, format, two_dimensions, max,
                       x_titles, y_titles, add_search_bar, contrasts = unique(dt$contrast)) {

  shared_df <- SharedData$new(dt, key = ~id, group = "A")
  filter_box <- filter_select("id_search", "Search by ID:", shared_df, ~id)
  showlegends <- c(TRUE, rep(FALSE, length(contrasts) - 1))
  names(showlegends) <- contrasts
  zerolines <- two_dimensions
  gg2 <- lapply(contrasts, function(cont) {
    message(cont)
    shared_df_sub <- SharedData$new(dt[contrast == cont,], key = ~id, group = "A")

    plot_ly(shared_df_sub, x = ~x_axis, y = ~y_axis, color = ~ Regulation, colors = color.values,
            text = ~trace, hoverinfo = "text", legendgroup = ~Regulation,
            type = "scattergl", showlegend = showlegends[cont],
            mode = "markers", marker = list(opacity = 0.5), name = ~ Regulation) %>%
      layout(xaxis = list(zerolinecolor = "red", zeroline = zerolines),
             yaxis = list(title = list(text = y_titles,
                                       font = list(color = "black", weight = "bold", size = 16)),
                          zerolinecolor = "red"),
             shapes = if (two_dimensions) {
               list(
                 list(type = "line", xref = "x", yref = "y",
                      x0 = max[1], x1 = max[2],
                      y0 = max[1], y1 = max[2],
                      line = list(color = "rgba(128, 128, 128, 0.5)", dash = "dash", width = 1)
                 )
               )}
      )
  })
  combined <- subplot(gg2, nrows = 1, shareX = TRUE, shareY = TRUE, titleX = FALSE, titleY = TRUE,
                      which_layout = "merge") %>%
    layout(
      annotations = list(
        list(x = 0.5, y = -0.15, text = x_titles[1], showarrow = FALSE, xref = "paper", yref = "paper",
             font = list(color = "black", weight = "bold", size = 16))),
      margin = list(b = 55),
      xaxis = list(range = c(max[1], max[2]))
    ) %>% plotly::config(toImageButtonOptions = list(format = format, filename = "RC_DEG_analysis"),
                         displaylogo = FALSE)

  widths <- cumsum(rep(1 / length(gg2), length(gg2)))
  center <- widths[1] / 2
  combined <- layout(combined, annotations = lapply(seq_along(contrasts), function(i) {
    list(x = widths[i] - center, y = 1.05, text = sub("Comparison: ", "", contrasts[i]),
         showarrow = FALSE, xref = 'paper', yref = 'paper', xanchor = 'center')
  }))
  attr(combined, "filter_box") <- filter_box
  if (add_search_bar) {
    combined <- DEG_add_search_bar(combined)
  }
  return(combined)
}

DEG_add_search_bar <- function(combined) {
  bscols(list(attr(combined, "filter_box"), tags$br(), combined), widths = 12)
}


DEG_plot_input_validation <- function() {
  with(rlang::caller_env(), {
    # colnames(dt)[colnames(dt) == "rna.lfc"] <- "rna"
    # colnames(dt)[colnames(dt) == "rfp.lfc"] <- "rfp"
    # colnames(dt)[colnames(dt) == "te.lfc"] <- "te"

    stopifnot(is(dt, "data.table"))
    if (is.character(xlim)) stopifnot(xlim %in% c("bidir.max", "auto"))
    if (is.character(ylim)) stopifnot(ylim %in% c("bidir.max", "auto"))
    # Remove this line in next bioc
    if("variable" %in% colnames(dt))
      colnames(dt) <- gsub("variable", "contrast", colnames(dt))
    columns_must_exists <- if (two_dimensions) {
      c("contrast", "Regulation", "id", "rna.lfc", "rfp.lfc", "te.lfc")
    } else c("contrast", "Regulation", "id", "LFC", "meanCounts")

    dt[, Regulation := factor(Regulation, levels = names(color.values))]

    if (!(all(columns_must_exists %in% colnames(dt))))
      stop("dt must minimally contain the columns: ",
           paste(columns_must_exists, collapse = ", "))
    if (!draw_non_regulated) dt <- dt[Regulation != "No change",]
    if (nrow(dt) == 0) {
      warning("dt input had no valid rows to plot")
    }
  })
}
