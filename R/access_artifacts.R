rc_access_path <- function(path, roots, context = rc_access_context()) {
  if (is.null(context)) return(invisible(TRUE))
  rc_access_check(context)
  path <- normalizePath(path, mustWork = TRUE)
  roots <- vapply(roots, normalizePath, character(1), mustWork = FALSE)
  if (any(!vapply(path, function(p) any(startsWith(p, paste0(roots, "/"))), logical(1))))
    stop("File unavailable or access denied.", call. = FALSE)
  invisible(TRUE)
}

rc_access_collection_columns <- function(table, context = rc_access_context()) {
  if (is.null(context)) return(table)
  rc_access_check(context)
  if ("library" %in% names(table)) return(table[library %in% context$runs])
  columns <- intersect(colnames(table), context$runs)
  if (!length(columns)) stop("Collection has no authorized libraries.", call. = FALSE)
  subset_collection_columns(table, columns)
}

rc_access_umap <- function(table, df, context = rc_access_context()) {
  if (is.null(context)) return(table)
  rc_access_check(context)
  column <- intersect(c("sample", "Run"), names(table))[1]
  if (is.na(column)) stop("UMAP lacks library identifiers; access denied.", call. = FALSE)
  table[as.character(table[[column]]) %in% ORFik::runIDs(df)]
}

rc_access_reference <- function(df, context = rc_access_context()) {
  if (is.null(context)) return(invisible(TRUE))
  rc_access_check(context)
  id <- context$catalog$id[context$catalog$experiment == ORFik::name(df)]
  allowed <- DBI::dbGetQuery(context$con,
    "SELECT dataset FROM reference_permissions WHERE dataset=? AND annotations=1", params = list(id))
  if (nrow(allowed) != 1L)
    stop("Shared reference annotations are not enabled for this dataset. Contact its administrator.", call. = FALSE)
  invisible(TRUE)
}
