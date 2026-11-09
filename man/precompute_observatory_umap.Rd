\name{precompute_observatory_umap}
\alias{precompute_observatory_umap}
\title{Precompute an Observatory default Plotly object}
\usage{precompute_observatory_umap(experiment, path = rc_umap_cache_file(experiment))}
\arguments{\item{experiment}{An ORFik collection experiment.}\item{path}{Output path; defaults to the collection reference UMAP directory.}}
\value{The output path, invisibly.}
\description{Build a cache for public sessions using the original 2D UMAP and tissue and cell-line coloring. The source checksum, collection Run order and Plotly version are verified before reuse.}
