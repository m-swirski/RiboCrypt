create_observatory_module <- function(organisms_df, samples_df) {}

# A vector with label of all available organisms
available_organisms <- c("")

# A vector with labels of all available samples
available_samples <- c("")

annotation_type_list <- c("")
library_type_list <- c("")
reading_frames <- c(1, 2, 3)

# A function that returns a list of genes for a given organism
# returns a vector of character vectors
get_gene_list <- function(organism) {}

# A function that return a list of isoforms for a given gene
# returns a vector of character vectors
get_isofrom_list <- function(gene) {}

# A function that returns position/positions of an isoform in the genome
# returns a GRangesList
get_region_for_isoform <- function(isoform) {}

extend_region_leader <- function(region) {}

extend_region_trailer <- function(region) {}

# A function that returns a list of samples available for a given organism
# returns a vector of character vectors
get_sample_list <- function(organism, library_types = library_type_list) {}

# TODO
# A function that returns a UMAP projection over samples for a given organisms
# returns a dataframe(TODO expand on the contents of dataframe)
# get_umap_data(oraganism, ?, ?, ?)

# TODO
# A function that returns coverage over a given region for selected samples
# returns a dataframe(TODO expand on the contents of the dataframe)
get_coverage_for_samples <- function(organism, samples, region) {}

# TODO
# It will possibly be a few different functions
# different for each normalization method
# normalize_coverage <- function(coverage) {}

# A function that smoothens coverage dataframe
aggregate_kmers <- function(coverage, kmer_size, aggregation_method) {}

# A function that adds reading frames metadata to coverage dataframe
mark_reading_frames <- function(coverage) {}

# A function that filters coverage dataframe, leaving only selected frames
subset_reading_frames <- function(
  coverage_with_frames, selected_frames = reading_frames
) {}

# A function that adds color metadata to coverage dataframe with marked reading frames
color_reading_frames <- function(coverage_with_frames) {}

# A function that smoothens coverage dataframe with marked reading frames
aggregate_kmers_with_reading_frames <- function(
  coverage_with_frames, kmer_size, aggregation_method
) {}

# TODO
# A function that returns annotations
# (Genes/Transcripts, Translons, uORFs, PhyloP, Mappability)
# returns a dataframe(TODO expand on the contents of the dataframe)
get_annotations <- function(
  organism, region, annotation_types = annotation_type_list
) {}

# TODO
# Figure out what dependencies a function like that should have
# When should it happen? Should it change the region before getting coverage?
# Or maybe it should happen as a post processing step on region/annotations?
# collapse_introns <- function() {}

create_coverage_track <- function(
  coverage,
  name,
  scale_transform = identity
) {}

create_annotation_track <- function(annotations) {}