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
  } else stop("Invalid Differential method selected")
  return(append_gene_symbols(dt, symbols_dt))
}

append_gene_symbols <- function(dt, symbols_dt) {
  stopifnot(is(dt, "data.table"))
  stopifnot(!is.null(dt$id))

  dt_with_symbols <- copy(dt)
  if (length(symbols_dt) > 0 & nrow(symbols_dt) > 0) {
    tx_column <- grep("tx|transcript|value", colnames(symbols_dt), ignore.case = TRUE, value = TRUE)
    if (length(tx_column) != 1) {
      if (length(tx_column) == 0) {
        warning("Could not find any column of symbols table with tx, transcript or value in column name, please rename/add")
      } else{
        warning("Found multiple columns of symbols table with tx, transcript or value in column name, please rename/remove")
      }
    }

    dt_with_symbols <- data.table::merge.data.table(dt, symbols_dt, by.x = "id", by.y = tx_column, sort = FALSE)
    if (!is.null(dt_with_symbols$external_gene_name)) {
      dt_with_symbols[, id := paste0(id, "(", external_gene_name, ")")]
    } else if (!is.null(dt_with_symbols$label)) {
      dt_with_symbols[, id := paste0(id, "(", sub("-.*", "", label), ")")]
    } else warning("Could not find any column of symbols table with external_gene_name or label in column name, please rename/add")
  }
  return(dt_with_symbols)
}

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
DEG_plot <- function(dt, draw_non_regulated = FALSE,
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
  contrasts <- unique(dt$contrast)
  id <- NULL
  dt <- highlight_key(dt, ~id, "Select a transcript")

  if (two_dimensions) {
    gg <- ggplot(dt, aes(x = rna, y = rfp, color = Regulation,
                         frame = contrast, id = id)) +
      geom_vline(aes(xintercept = 0), color = "red") +
      geom_abline(slope = 1, intercept = 0, color = "grey", linetype = 2)
  } else { # One dimensional
    gg <- ggplot(dt, aes(x = log2(meanCounts), y = LFC, color = Regulation,
                         frame = contrast, id = id))
    if (length(contrasts) == 1) gg <- gg + ggtitle(contrasts)
  }
   gg <- gg +
    geom_hline(aes(yintercept = 0), color = "red") +
    geom_point() + xlab(xlab) + ylab(ylab) +
    scale_color_manual(values = color.values)
   if (!is.null(xlim)) gg <- gg + xlim(xlim)
   if (!is.null(ylim)) gg <- gg + ylim(ylim)

   hovertip <- "id"

   #vis <- ifelse(unique(dt$Regulation) == "", TRUE, "legendonly") #add_trace(visible = vis)
  #removed partial_bundle() - it's incompatible with shiny
   select <- highlight(
     ggplotly(gg, tooltip = hovertip) %>%
       layout(autosize = TRUE) %>% plotly::config(
         toImageButtonOptions = list(
           format = format),
         displaylogo = FALSE),
     on = c('plotly_selected'), off = c('plotly_deselect'),
     selectize = TRUE, persistent = FALSE
   )
   if (length(contrasts) == 1) {
     select <- select %>% toWebGL()
   }

   return(select)
}

DEG_plot_input_validation <- function() {
  with(rlang::caller_env(), {
    colnames(dt)[colnames(dt) == "rna.lfc"] <- "rna"
    colnames(dt)[colnames(dt) == "rfp.lfc"] <- "rfp"
    colnames(dt)[colnames(dt) == "te.lfc"] <- "te"

    stopifnot(is(dt, "data.table"))
    if (is.character(xlim)) stopifnot(xlim %in% c("bidir.max", "auto"))
    if (is.character(ylim)) stopifnot(ylim %in% c("bidir.max", "auto"))
    # Remove this line in next bioc
    if("variable" %in% colnames(dt))
      colnames(dt) <- gsub("variable", "contrast", colnames(dt))
    columns_must_exists <- if (two_dimensions) {
      c("contrast", "Regulation", "id", "rna", "rfp", "te")
    } else c("contrast", "Regulation", "id", "LFC", "meanCounts")

    if (!(all(columns_must_exists %in% colnames(dt))))
      stop("dt must minimally contain the columns: ",
           paste(columns_must_exists, collapse = ", "))
    if (!draw_non_regulated) dt <- dt[Regulation != "No change",]
    if (nrow(dt) == 0) {
      warning("dt input had no valid rows to plot")
    }
  })
}
