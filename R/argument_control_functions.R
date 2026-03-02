# Controller for input validation
multiOmicsController <- function() {
  with(rlang::caller_env(), {
    annotation_list <- annotation_controller(
      df, display_range,
      annotation, annotation_names,
      leader_extension, trailer_extension, viewMode
    )
    display_range <- annotation_list$display_range
    annotation <- annotation_list$annotation
  })
}

# Controller for view input validation
multiOmicsControllerView <- function(
  use_fst = FALSE, selected_libraries = list()
) {
  with(rlang::caller_env(), {
    if (!is(reads, "list")) reads <- list(reads)

    if (length(withFrames) == 0) withFrames <- FALSE
    if (!(length(withFrames) %in% c(1, length(reads), length(selected_libraries)))) stop("length of withFrames must be 0, 1 or the same as reads list")
    if (length(withFrames) == 1) withFrames <- rep(withFrames, length(reads))

    if (length(colors) == 0 && !use_fst) colors <- 1:length(reads)
    if (length(colors) == 0 && use_fst) colors <- length(selected_libraries)

    if (!(length(colors) %in% c(1, length(reads), length(selected_libraries)))) stop("length of colors must be 0, 1 or the same as reads list")

    if (length(colors) == 1 && !use_fst) colors <- rep(colors, length(reads))
    if (length(colors) == 1 && use_fst) colors <- rep(colors, length(selected_libraries))

    if (length(kmers) == 0) kmers <- 1
    if (!(length(kmers) %in% c(1, length(reads), length(selected_libraries)))) stop("length of kmers must be 0, 1 or the same as reads list")

    if (length(kmers) == 1 && !use_fst) kmers <- rep(kmers, length(reads))
    if (length(kmers) == 1 && use_fst) kmers <- rep(kmers, length(selected_libraries))

    if (length(ylabels) == 0 && !use_fst) ylabels <- as.character(1:length(reads))
    if (length(ylabels) == 0 && use_fst) ylabels <- as.character(1:length(selected_libraries))

    if (!(length(ylabels) %in% c(1, length(reads), length(selected_libraries)))) stop("length of ylabels must be 0, 1 or the same as reads list")

    if (length(ylabels) == 1 && !use_fst) ylabels <- rep(ylabels, length(reads))
    if (length(ylabels) == 1 && use_fst) ylabels <- rep(ylabels, length(selected_libraries))

    ylabels_full_name <- ylabels

    if (length(ylabels) > 5 && !use_fst) ylabels <- 1:length(reads)
    if (length(ylabels) > 5 && use_fst) ylabels <- 1:length(selected_libraries)

    if (length(lib_proportions) == 0) lib_proportions <- 1
    if (!(length(lib_proportions) %in% c(1, length(reads), length(selected_libraries)))) stop("length of lib_proportions must be 0, 1 or the same as reads list")
    if (length(lib_proportions) == 1) {
      lib_proportions <- if (frames_type == "animate") {
        lib_proportions
      } else {
        if (!use_fst) rep(lib_proportions, length(reads)) else rep(lib_proportions, length(selected_libraries))
      }
    }
    lib_proportions <- lib_proportions / sum(lib_proportions)
    if (length(annotation_proportions) == 0) {
      if (display_sequence %in% c("none", FALSE)) {
        annotation_proportions <- c(0.35, 0.65)
      } else if (viewMode == "genomic") {
        lib_to_annotation_proportions <- c(0.6, 0.4)
        annotation_proportions <- c(0.2, 0.5, 0.3)
      } else { # tx coordinates
        annotation_proportions <- c(0.2, 0.2, 0.6)
        if (bottom_panel$annotation_layers > 1) {
          annotation_proportions[2] <- annotation_proportions[2] * sqrt(bottom_panel$annotation_layers)
        }
      }
    } else {
      annotation_proportions <- annotation_proportions / sum(annotation_proportions)
    }
    proportions <- c(lib_proportions * lib_to_annotation_proportions[1], annotation_proportions * lib_to_annotation_proportions[2])
    if (summary_track) proportions <- c(0.2, proportions * 0.8)

    if (bottom_panel$ncustom > 0) {
      proportions <- c(proportions, rep(0.12, bottom_panel$ncustom))
    }
    proportions <- proportions / sum(proportions)
  })
}

# Controller for correct annotation definition
annotation_controller <- function(df, display_range, annotation, annotation_names = NULL,
                                  leader_extension, trailer_extension, viewMode) {
  # Load ranges if defined as character
  if (is(display_range, "character")) {
    display_range <- loadRegion(df, part = "tx", names.keep = display_range)
    if (length(display_range) == 0) stop("display_range specified as name of non existing transcript,
                                       did you specify a name from a different genome?")
  }
  if (is(annotation, "character")) {
    annotation <- loadRegion(df, part = annotation)
    if (length(annotation) == 0) stop("annotation specified on valid transcript,
                                    but does not contain this region, did you specify cds on
                                    a non coding RNA or leader for mRNA without defined leader?")
  }

  seqlevels(display_range) <- seqlevels(annotation)
  display_range <- GRangesList(display_range)

  if (!is.null(leader_extension) && is.numeric(leader_extension) && leader_extension != 0) {
    display_range <- extendLeaders(display_range, leader_extension)
  }
  if (!is.null(trailer_extension) && is.numeric(trailer_extension) && trailer_extension != 0) {
    display_range <- extendTrailers(display_range, trailer_extension)
  }

  if (!is.null(annotation_names)) {
    if (length(annotation_names) == 1) {
      if (annotation_names %in% annotation) {
        names(annotation) <- mcols(annotation)[[annotation_names]]
      } else {
        stop(wmsg("wrong annotation_names argument"))
      }
    } else if (length(annotation_names) == length(annotation)) {
      names(annotation) <- annotation_names
    } else {
      stop(wmsg("wrong annotation_names argument"))
    }
  }
  return(list(display_range = display_range, annotation = annotation))
}
