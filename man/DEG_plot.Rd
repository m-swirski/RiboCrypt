% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DEG.R
\name{DEG_plot}
\alias{DEG_plot}
\title{Differential expression plots (1D or 2D)}
\usage{
DEG_plot(
  dt,
  draw_non_regulated = TRUE,
  add_search_bar = TRUE,
  xlim = ifelse(two_dimensions, "bidir.max", "auto"),
  ylim = "bidir.max",
  xlab = ifelse(two_dimensions, "RNA fold change (log2)", "Mean counts (log2)"),
  ylab = ifelse(two_dimensions, "RFP fold change (log2)", "Fold change (log2)"),
  two_dimensions = ifelse("LFC" \%in\% colnames(dt), FALSE, TRUE),
  color.values = c(`No change` = "black", Significant = "red", Buffering = "purple",
    `mRNA abundance` = "darkgreen", Expression = "blue", Forwarded = "yellow", Inverse =
    "aquamarine", Translation = "orange4"),
  format = "png"
)
}
\arguments{
\item{dt}{a data.table with results from a differential
expression run. Normally from: \code{ORFik::DTEG.analysis(df1, df2)}}

\item{draw_non_regulated}{logical, default TRUE
Should non-regulated rows be included in the plot?
Will make the plot faster to render if skipped (FALSE)}

\item{add_search_bar}{logical, default TRUE. Add a crosstalk search bar
to search for genes in the plot}

\item{xlim}{numeric vector or character preset, default:
\code{ifelse(two_dimensions, "bidir.max", "auto")}
(Equal in both + / - direction, using max value + 0.5 of
meanCounts(in 1d) / rna(in 2d) column of dt).
If you want ggplot to decide limit, set to "auto". For numeric vector,
specify min and max x limit: like c(-5, 5)}

\item{ylim}{numeric vector or character preset, default: "bidir.max"
(Equal in both + / - direction, using max value + 0.5 of
LFC(in 1d) / rfp(in 2d) column of dt).
If you want ggplot to decide limit, set to "auto". For numeric vector,
specify min and max x limit: like c(-5, 5)}

\item{xlab}{character, default:
\code{ifelse(two_dimensions, "RNA fold change (log2)", "Mean counts (log2)")}}

\item{ylab}{character, default:
\code{ifelse(two_dimensions, "RFP fold change (log2)",  "Fold change (log2)")}}

\item{two_dimensions}{logical, default:
\code{ifelse("LFC" \%in\% colnames(dt), FALSE, TRUE)}
Is this two dimensional, like: Ribo-seq vs RNA-seq. Alternative, FALSE: Then
Log fold change vs mean counts}

\item{color.values}{named character vector, default: \code{c(
"No change" = "black", "Significant" = "red",
"Buffering" = "purple", "mRNA abundance" = "darkgreen",
"Expression" = "blue", "Forwarded" = "yellow",
"Inverse" = "aquamarine", "Translation" = "orange4")}}

\item{format}{character, default "png". Format for plotly bar.}
}
\value{
plotly object or crosstalk bscols if add_search_bar is TRUE.
}
\description{
Gives you interactive 1D or 2D DE plots
}
\examples{
# Load experiment
df <- ORFik.template.experiment()
df_rna <- df[df$libtype == "RNA",]
# 1 Dimensional analysis
dt <- DEG.analysis(df_rna)
dt$Regulation[1] <- "Significant" # Fake sig level
DEG_plot(dt, draw_non_regulated = TRUE)
# 2 Dimensional analysis
df_rfp <- df[df$libtype == "RFP",]
dt_2d <- DTEG.analysis(df_rfp, df_rna, output.dir = NULL)
dt_2d$Regulation[4] <- "Translation" # Fake sig level
dt_2d$rfp.lfc[4] <- -0.3 # Fake sig level
dt_2d$Regulation[5] <- "Buffering" # Fake sig level
dt_2d$rna.lfc[5] <- -0.3 # Fake sig level
DEG_plot(dt_2d, draw_non_regulated = TRUE)
}
