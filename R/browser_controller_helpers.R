#' @noRd
sanitize_library_selection_group <- function(selection) {
  selection <- as.character(selection)
  unique(selection[!is.na(selection) & nzchar(selection)])
}

#' @noRd
normalize_library_selections <- function(library_selections) {
  if (is.null(library_selections)) return(list())

  lib_sel <- lapply(library_selections, sanitize_library_selection_group)
  lib_sel[lengths(lib_sel) > 0]
}

#' @noRd
match_library_selections_to_runs <- function(library_selections, available_runs) {
  lib_sel <- lapply(library_selections, intersect, y = available_runs)
  lib_sel <- lib_sel[lengths(lib_sel) > 0]

  list(
    selections = lib_sel,
    selection_indices = lapply(lib_sel, match, table = available_runs)
  )
}

#' @noRd
selection_display_label <- function(selection_id, library_selection_labels = NULL) {
  label <- library_selection_labels[[selection_id]]
  if (is.null(label) || label == "" || label == selection_id) {
    selection_id
  } else {
    paste(selection_id, label, sep = " - ")
  }
}

#' @noRd
apply_library_selection_labels <- function(library_selections, library_selection_labels = NULL) {
  if (is.null(library_selection_labels)) return(library_selections)

  selection_ids <- names(library_selections)
  names(library_selections) <- vapply(
    selection_ids,
    selection_display_label,
    character(1),
    library_selection_labels = library_selection_labels
  )
  library_selections
}

#' @noRd
rename_all_merged_library_selections <- function(library_selections, collection_path) {
  fst_index <- collection_path
  file_forward <- fst::read_fst(collection_path)[1, ]$file_forward
  file_forward <- file.path(dirname(fst_index), basename(file_forward))
  megafst_samples <- length(fst::metadata_fst(file_forward)$columnNames)
  group_is_all <- lengths(library_selections) == megafst_samples
  if (any(group_is_all)) names(library_selections)[group_is_all] <- "All merged"

  message("Number of runs used: ", length(unlist(library_selections)),
          " (", paste(lengths(library_selections), collapse = ", "), ")")

  library_selections
}

#' @noRd
collection_profile_template <- function(reads_mat) {
  data.table::data.table(
    position = seq_len(nrow(reads_mat)),
    library = as.factor(rep_len(1, length.out = nrow(reads_mat))),
    frame = as.factor(rep_len(1:3, length.out = nrow(reads_mat)))
  )
}

#' @noRd
build_collection_profile <- function(selection_idx, reads_mat, profile_template, kmer) {
  cbind(
    data.table::data.table(
      count = matrixStats::rowSums2(
        reads_mat,
        cols = selection_idx,
        na.rm = TRUE
      )
    ),
    profile_template
  ) |>
    smoothenMultiSampCoverage(
      kmer,
      kmers_type = "mean",
      split_by_frame = TRUE
    )
}

#' @noRd
build_collection_profiles <- function(selection_indices, reads_mat, kmer) {
  profile_template <- collection_profile_template(reads_mat)
  lapply(
    selection_indices,
    build_collection_profile,
    reads_mat = reads_mat,
    profile_template = profile_template,
    kmer = kmer
  )
}

#' @noRd
collection_library_hash_suffix <- function(library_selections) {
  paste(
    vapply(names(library_selections), function(selection_id) {
      paste(selection_id, paste(library_selections[[selection_id]], collapse = ","), sep = ":")
    }, character(1)),
    collapse = "|group|"
  )
}

#' @noRd
build_single_browser_track_plot <- function(profile, with_frames, frame_color, color,
                                            ylabel, ylabel_full_name, lines,
                                            frames_type, lib_index, total_libs,
                                            templates = NULL) {
  createSinglePlot(
    profile,
    with_frames,
    frame_color,
    color,
    ylabel,
    ylabel_full_name,
    lines,
    type = frames_type,
    lib_index,
    total_libs,
    templates = templates
  )
}

#' @noRd
build_browser_track_plots <- function(profiles, withFrames, frame_colors, colors,
                                      ylabels, ylabels_full_name, lines,
                                      frames_type, BPPARAM, templates = NULL) {
  total_libs <- length(profiles)
  plot_args <- list(
    profiles,
    withFrames,
    frame_colors,
    colors,
    ylabels,
    ylabels_full_name,
    seq_len(total_libs)
  )

  if (is(BPPARAM, "SerialParam")) {
    return(do.call(
      mapply,
      c(
        list(
          FUN = build_single_browser_track_plot
        ),
        plot_args,
        list(
          MoreArgs = list(
            lines = lines,
            frames_type = frames_type,
            total_libs = total_libs,
            templates = templates
          ),
          SIMPLIFY = FALSE
        )
      )
    ))
  }

  do.call(
    bpmapply,
    c(
      list(FUN = build_single_browser_track_plot),
      plot_args,
      list(
        MoreArgs = list(
          lines = lines,
          frames_type = frames_type,
          total_libs = total_libs,
          templates = templates
        ),
        SIMPLIFY = FALSE,
        BPPARAM = BPPARAM
      )
    )
  )
}

#' @noRd
build_single_browser_profile <- function(read, with_frames, kmer, log_scale,
                                         display_range, kmers_type, frames_type,
                                         frames_subset) {
  getProfileWrapper(
    display_range,
    read,
    with_frames,
    kmer,
    log_scale,
    kmers_type,
    type = frames_type,
    frames_subset = frames_subset
  )
}

#' @noRd
build_browser_profiles <- function(display_range, reads, kmers, kmers_type,
                                   frames_type, frames_subset, withFrames,
                                   log_scale, BPPARAM) {
  profile_args <- list(reads, withFrames, kmers, log_scale)

  if (is(BPPARAM, "SerialParam")) {
    return(do.call(
      mapply,
      c(
        list(FUN = build_single_browser_profile),
        profile_args,
        list(
          MoreArgs = list(
            display_range = display_range,
            kmers_type = kmers_type,
            frames_type = frames_type,
            frames_subset = frames_subset
          ),
          SIMPLIFY = FALSE
        )
      )
    ))
  }

  do.call(
    bpmapply,
    c(
      list(FUN = build_single_browser_profile),
      profile_args,
      list(
        MoreArgs = list(
          display_range = display_range,
          kmers_type = kmers_type,
          frames_type = frames_type,
          frames_subset = frames_subset
        ),
        SIMPLIFY = FALSE,
        BPPARAM = BPPARAM
      )
    )
  )
}
