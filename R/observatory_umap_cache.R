# Cached Plotly objects preserve Shiny's onRender selection hook, which is
# attached separately for each session after loading. Never use a shared cache
# in authenticated sessions: the current collection may be permission-filtered.
rc_umap_cache_file <- function(df) {
  staging <- getOption("ribocrypt.umap_cache_dir",NULL)
  if (!is.null(staging)) return(file.path(staging,paste0(ORFik::name(df),"_observatory_default.rds")))
  safe_name <- gsub("[^A-Za-z0-9_.-]","_",ORFik::name(df))
  file.path(ORFik::refFolder(df),"UMAP",paste0("observatory_plotly_",safe_name,"_tissue_cellline_6d_v2.rds"))
}
rc_umap_source_file <- function(df) file.path(ORFik::refFolder(df),"UMAP","UMAP_by_gene_counts_6d.fst")
rc_umap_render_key <- function(df, embedding = "2d") {
  dir <- file.path(ORFik::refFolder(df),"UMAP")
  files <- c(file.path(dir,"UMAP_by_gene_counts.fst"))
  if (identical(embedding,"6d")) files <- c(files,
    file.path(dir,"UMAP_by_gene_counts_6d.fst"),
    file.path(dir,"UMAP_by_gene_counts_inputs.rds"),
    file.path(dir,"UMAP_by_gene_counts_6d_inputs.rds"))
  annotation <- file.path(dir, paste0("UMAP_", embedding, "_axis_annotations.csv"))
  if (identical(ORFik::name(df), "all_samples-Homo_sapiens")) files <- c(files, annotation)
  hashes <- vapply(files,function(x) if (file.exists(x)) unname(tools::md5sum(x)) else NA_character_,character(1))
  # Include the session's actual Run order: access-filtered collections can
  # share a name while exposing different libraries to different users.
  list(hashes=hashes, mtimes = vapply(files, function(x) if (file.exists(x)) as.numeric(file.info(x)$mtime) else NA_real_, numeric(1)), runs=as.character(ORFik::runIDs(df)))
}
rc_umap_cache_source_hash <- function(df) {
  dir <- file.path(ORFik::refFolder(df), "UMAP")
  files <- file.path(dir, c("UMAP_by_gene_counts_6d.fst", "UMAP_by_gene_counts.fst",
                            "UMAP_by_gene_counts_6d_inputs.rds", "UMAP_by_gene_counts_inputs.rds"))
  if (any(!file.exists(files))) return(NA_character_)
  paste(unname(tools::md5sum(files)), collapse = ":")
}
rc_umap_cache_valid <- function(x, source_md5, runs, plotly_version) {
  !is.na(source_md5) && is.list(x) && identical(x$schema,2L) && identical(x$source_md5,source_md5) &&
    identical(x$plotly_version,plotly_version) && identical(x$runs,runs) && inherits(x$plot,"plotly")
}

rc_umap_cached_plot <- function(experiment, color_by, embedding, axes) {
  if (is.null(experiment) || !is.null(rc_access_context()) || !isTRUE(getOption("ribocrypt.umap_cache",TRUE)) ||
      !identical(as.character(color_by),c("tissue","cell_line")) || !identical(embedding,"6d") ||
      !identical(as.integer(axes),c(1L,2L))) return(NULL)
  path <- rc_umap_cache_file(experiment); source <- rc_umap_source_file(experiment)
  if (!file.exists(path) || !file.exists(source)) return(NULL)
  tryCatch({
    x <- readRDS(path)
    if (!rc_umap_cache_valid(x, rc_umap_cache_source_hash(experiment), as.character(ORFik::runIDs(experiment)),
                             as.character(utils::packageVersion("plotly")))) return(NULL)
    x$plot
  }, error=function(e) NULL)
}

#' Precompute an Observatory default Plotly object for one collection
#'
#' The cache is used only for public sessions, the six-dimensional projection and
#' the default tissue/cell-line coloring. Source checksum, collection Run order,
#' and Plotly version are checked when it is loaded. Stale caches are ignored.
#' @param experiment An ORFik collection experiment.
#' @param path Output path; defaults to the collection reference UMAP directory.
#' @return The output path, invisibly.
#' @export
precompute_observatory_umap <- function(experiment, path=rc_umap_cache_file(experiment)) {
  if (!is.null(rc_access_context())) stop("Shared UMAP caches cannot be built within an access-filtered session")
  source <- rc_umap_source_file(experiment)
  if (!file.exists(source)) stop("No six-dimensional UMAP available: ",source)
  before <- rc_umap_cache_source_hash(experiment)
  if (is.na(before)) stop("UMAP cache requires both coordinate tables and input provenance")
  d <- load_data_umap_internal(ORFik::name(experiment), color.by=c("tissue", "cell_line"), embedding="6d")
  runs <- as.character(ORFik::runIDs(experiment))
  if (!identical(as.character(d$sample),runs)) stop("UMAP sample order differs from the collection; rebuild the UMAP first")
  d[, color_column := paste(tissue,cell_line,sep=" | ")]
  p <- umap_plot(d,color.by=c("tissue","cell_line"),template=observatoryUmapPlotlyTemplate())
  if (!identical(before,rc_umap_cache_source_hash(experiment))) stop("UMAP source changed during cache creation")
  x <- list(schema=2L, source_md5=before, plotly_version=as.character(utils::packageVersion("plotly")), runs=runs, plot=p)
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  tmp <- tempfile(pattern=".umap-cache-",tmpdir=dirname(path)); on.exit(unlink(tmp),add=TRUE)
  saveRDS(x,tmp,compress=FALSE)
  if (!file.rename(tmp,path)) stop("Could not publish the precomputed UMAP cache")
  invisible(path)
}

# Interpretation is tied to the exact human 2D coordinates, never inferred
# from a different embedding or transferred to another organism.
rc_umap_cached_axis_titles <- function(plot, experiment, source) {
  if (!identical(ORFik::name(experiment), "all_samples-Homo_sapiens")) return(plot)
  annotation <- file.path(dirname(source), "UMAP_2d_axis_annotations.csv")
  if (!file.exists(annotation) || !file.exists(source)) return(plot)
  tryCatch({
    if (file.info(annotation)$mtime < file.info(source)$mtime) return(plot)
    a <- utils::read.csv(annotation, stringsAsFactors = FALSE, check.names = FALSE)
    if (!all(c("axis", "label", "source_md5") %in% names(a)) ||
        nrow(a) != 2L || anyDuplicated(a$axis) ||
        !setequal(a$axis, 1:2) || anyNA(a$label) || any(!nzchar(a$label)) ||
        !all(a$source_md5 == unname(tools::md5sum(source)))) return(plot)
    a <- a[match(1:2, a$axis), , drop = FALSE]
    plot$x$layout$xaxis$title <- paste0("UMAP 1 (Enriched: ", a$label[1], ")")
    plot$x$layout$yaxis$title <- paste0("UMAP 2 (Enriched: ", a$label[2], ")")
    plot
  }, error = function(e) plot)
}

# 6D is a separate embedding; annotations are versioned against its own file.
rc_umap_6d_axis_titles <- function(plot, experiment, axes) {
  if (is.null(experiment) || !identical(ORFik::name(experiment), "all_samples-Homo_sapiens")) return(plot)
  source <- file.path(ORFik::refFolder(experiment), "UMAP", "UMAP_by_gene_counts_6d.fst")
  annotation <- file.path(dirname(source), "UMAP_6d_axis_annotations.csv")
  if (!file.exists(annotation) || !file.exists(source)) return(plot)
  tryCatch({
    if (file.info(annotation)$mtime < file.info(source)$mtime) return(plot)
    a <- utils::read.csv(annotation, stringsAsFactors = FALSE, check.names = FALSE)
    if (!all(c("axis", "label", "source_md5") %in% names(a)) ||
        nrow(a) != 6L || anyDuplicated(a$axis) || !setequal(a$axis, 1:6) ||
        anyNA(a$label) || any(!nzchar(a$label)) ||
        !all(a$source_md5 == unname(tools::md5sum(source)))) return(plot)
    labels <- a$label[match(as.integer(axes), a$axis)]
    if (length(labels) != 2L || anyNA(labels)) return(plot)
    # layout() handles both fully built templates and normal plotly objects.
    plotly::layout(plot, xaxis = list(title = paste0("UMAP ", axes[1], " (Enriched: ", labels[1], ")")),
                   yaxis = list(title = paste0("UMAP ", axes[2], " (Enriched: ", labels[2], ")")))
  }, error = function(e) plot)
}
