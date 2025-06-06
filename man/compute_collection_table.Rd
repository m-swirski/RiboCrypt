% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/collection_helpers.R
\name{compute_collection_table}
\alias{compute_collection_table}
\title{Get collection table normalized in wide format}
\usage{
compute_collection_table(
  path,
  lib_sizes,
  df,
  metadata_field,
  normalization,
  kmer,
  metadata,
  min_count = 0,
  format = "wide",
  value.var = "logscore",
  as_list = FALSE,
  subset = NULL,
  group_on_tx_tpm = NULL,
  split_by_frame = FALSE,
  ratio_interval = NULL,
  decreasing_order = FALSE
)
}
\arguments{
\item{path}{the path to gene counts}

\item{lib_sizes}{named integer vector, default NULL. If given will do a pre tpm normalization
for full library sizes}

\item{df}{the ORFik experiment to load the precomputed collection from.
It must also have defined runIDs() for all samples.}

\item{metadata_field}{the column name in metadata, to select to group on.}

\item{normalization}{a character string, which mode, for options see RiboCrypt:::normalizations}

\item{kmer}{integer, default 1L (off), if > 1 will smooth out signal with sliding window size kmer.}

\item{metadata}{a data.table of metadata, must contain the Run column to
select libraries.}

\item{min_count}{integer, default 0. Minimum counts of coverage over transcript
to be included.}

\item{format}{character, default "wide", alternative "long". The format of
the table output.}

\item{value.var}{which column to use as scores, default "logscore"}

\item{as_list}{logical, default FALSE. Return as list of size 2,
count data.table and metadata data.table Set to TRUE if you need metadata
subset (needed if you subset the table, to get correct matching)}

\item{subset}{numeric vector, positional interval to subset, must be <= size of
whole region.}

\item{group_on_tx_tpm}{numeric vector, default NULL.
tpm values per libraries. Either for that gene or some other gene.}

\item{split_by_frame}{logical, default FALSE
For kmer sliding window, should it split by frame}

\item{ratio_interval}{numeric vector of size 2 or 4, default NULL.
If 2, means you should
sort libraries on coverage in that region. If 4, means to sort on ratio
of that region in this gene vs the other region in another gene.}

\item{decreasing_order}{logical, default FALSE. Sort you ordering vector from lowest (default).
If TRUE, sort from highest downwards.}
}
\value{
a data.table in long or wide (default) format, if as list, it is a
list of size 2 (see argument as_list)
}
\description{
Get collection table normalized in wide format
}
