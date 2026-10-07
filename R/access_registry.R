# Optional administrator-owned registry; no password material belongs here.
access_scalar <- function(value, name) {
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !nzchar(trimws(value)) || grepl("[[:cntrl:]]", value))
    stop(name, " must be a nonempty string without control characters.", call. = FALSE)
  value
}

access_flag <- function(value, name) {
  if (!is.logical(value) || length(value) != 1L || is.na(value))
    stop(name, " must be TRUE or FALSE.", call. = FALSE)
  as.integer(value)
}

#' Open the Optional Account and Dataset Registry
#' @param path SQLite database path, or `:memory:` for tests.
#' @return An open DBI connection. The caller must disconnect it.
#' @export
ribocrypt_access_db <- function(path) {
  access_scalar(path, "path")
  if (!requireNamespace("DBI", quietly = TRUE) ||
      !requireNamespace("RSQLite", quietly = TRUE))
    stop("Account registry requires DBI and RSQLite.", call. = FALSE)
  new_file <- path != ":memory:" && !file.exists(path)
  con <- DBI::dbConnect(RSQLite::SQLite(), path)
  if (new_file) Sys.chmod(path, "0600")
  tryCatch(access_schema(con), error = function(e) { DBI::dbDisconnect(con); stop(e) })
  con
}

access_schema <- function(con) {
  DBI::dbExecute(con, "PRAGMA foreign_keys = ON")
  DBI::dbExecute(con, "PRAGMA busy_timeout = 5000")
  statements <- c(
    "CREATE TABLE IF NOT EXISTS accounts (issuer TEXT NOT NULL, subject TEXT NOT NULL, first_seen REAL NOT NULL, last_seen REAL NOT NULL, PRIMARY KEY (issuer,subject))",
    "CREATE TABLE IF NOT EXISTS workspaces (id TEXT PRIMARY KEY, label TEXT NOT NULL)",
    "CREATE TABLE IF NOT EXISTS members (issuer TEXT NOT NULL, subject TEXT NOT NULL, workspace TEXT NOT NULL REFERENCES workspaces(id), PRIMARY KEY (issuer, subject, workspace))",
    "CREATE TABLE IF NOT EXISTS datasets (id TEXT PRIMARY KEY, experiment TEXT NOT NULL, directory TEXT NOT NULL, public INTEGER NOT NULL CHECK(public IN (0,1)))",
    "CREATE TABLE IF NOT EXISTS grants (workspace TEXT NOT NULL REFERENCES workspaces(id), dataset TEXT NOT NULL REFERENCES datasets(id), download INTEGER NOT NULL CHECK(download IN (0,1)), PRIMARY KEY (workspace, dataset))",
    "CREATE TABLE IF NOT EXISTS reference_permissions (dataset TEXT PRIMARY KEY REFERENCES datasets(id), annotations INTEGER NOT NULL CHECK(annotations IN (0,1)))",
    "CREATE INDEX IF NOT EXISTS dataset_visibility ON datasets(public)",
    "CREATE INDEX IF NOT EXISTS grant_dataset ON grants(dataset)"
  )
  DBI::dbWithTransaction(con, for (sql in statements) DBI::dbExecute(con, sql))
  invisible(con)
}

#' Register a Workspace
#' @param con Registry connection.
#' @param id Immutable workspace ID.
#' @param label Display name.
#' @export
ribocrypt_access_workspace <- function(con, id, label = id) {
  params <- list(access_scalar(id, "id"), access_scalar(label, "label"))
  DBI::dbExecute(con, "INSERT INTO workspaces (id,label) VALUES (?,?) ON CONFLICT(id) DO UPDATE SET label=excluded.label", params = params)
  invisible(id)
}

#' Assign or Remove Workspace Membership
#' @param con Registry connection.
#' @param issuer OIDC issuer URL.
#' @param subject Immutable OIDC subject, not an email address.
#' @param workspace Workspace ID.
#' @param enabled Whether membership is enabled.
#' @export
ribocrypt_access_member <- function(con, issuer, subject, workspace, enabled = TRUE) {
  params <- list(access_scalar(issuer, "issuer"), access_scalar(subject, "subject"),
                 access_scalar(workspace, "workspace"))
  sql <- if (access_flag(enabled, "enabled"))
    "INSERT OR IGNORE INTO members VALUES (?,?,?)" else
    "DELETE FROM members WHERE issuer=? AND subject=? AND workspace=?"
  DBI::dbExecute(con, sql, params = params)
  invisible(workspace)
}

#' Register an Immutable Dataset
#' @param con Registry connection.
#' @param id Immutable dataset/version ID.
#' @param experiment ORFik experiment name (not a file path).
#' @param directory Existing experiment directory.
#' @param public Whether anonymous users may read and download this dataset.
#' @param reference_annotations Whether reference-wide translon predictions and
#'   protein structures are approved for all readers of this dataset.
#' @export
ribocrypt_access_dataset <- function(con, id, experiment, directory, public = FALSE,
                                     reference_annotations = FALSE) {
  access_scalar(experiment, "experiment")
  if (grepl("[/\\\\]", experiment) || experiment %in% c(".", ".."))
    stop("experiment must be a name, not a path.", call. = FALSE)
  access_scalar(directory, "directory")
  if (!dir.exists(directory)) stop("Experiment directory does not exist.", call. = FALSE)
  params <- list(access_scalar(id, "id"), experiment,
                 normalizePath(directory, mustWork = TRUE), access_flag(public, "public"))
  annotations <- access_flag(reference_annotations, "reference_annotations")
  DBI::dbWithTransaction(con, {
    DBI::dbExecute(con, "INSERT INTO datasets (id,experiment,directory,public) VALUES (?,?,?,?)", params = params)
    DBI::dbExecute(con, "INSERT INTO reference_permissions VALUES (?,?)", params = list(id, annotations))
  })
  invisible(id)
}

#' Assign or Revoke a Dataset Grant
#' @param con Registry connection.
#' @param workspace Workspace ID.
#' @param dataset Dataset ID.
#' @param download Whether export is allowed in addition to reading.
#' @param enabled Whether the grant is enabled.
#' @export
ribocrypt_access_grant <- function(con, workspace, dataset, download = FALSE, enabled = TRUE) {
  params <- list(access_scalar(workspace, "workspace"), access_scalar(dataset, "dataset"))
  can_download <- access_flag(download, "download")
  sql <- if (access_flag(enabled, "enabled"))
    "INSERT INTO grants VALUES (?,?,?) ON CONFLICT(workspace,dataset) DO UPDATE SET download=excluded.download" else
    "DELETE FROM grants WHERE workspace=? AND dataset=?"
  if (enabled) params <- c(params, list(can_download))
  DBI::dbExecute(con, sql, params = params)
  invisible(dataset)
}

#' List Authorized Datasets
#' @param con Registry connection.
#' @param identity NULL for anonymous, or a verified issuer/subject list.
#' @param capability Either `read` or `download`.
#' @return A data frame of authorized immutable dataset IDs and storage mappings.
#' @export
ribocrypt_access_catalog <- function(con, identity = NULL, capability = c("read", "download")) {
  capability <- match.arg(capability)
  if (is.null(identity)) return(DBI::dbGetQuery(con, "SELECT * FROM datasets WHERE public=1 ORDER BY id"))
  params <- list(access_scalar(identity$issuer, "issuer"), access_scalar(identity$subject, "subject"))
  condition <- if (capability == "download") " AND g.download=1" else ""
  sql <- paste0("SELECT d.* FROM datasets d WHERE d.public=1 OR d.id IN (SELECT g.dataset FROM grants g JOIN members m ON m.workspace=g.workspace WHERE m.issuer=? AND m.subject=?", condition, ") ORDER BY d.id")
  DBI::dbGetQuery(con, sql, params = params)
}

#' Resolve a Dataset After Checking Current Permissions
#' @param con Registry connection.
#' @param dataset Immutable dataset ID, never an experiment name or path.
#' @param identity Verified identity, or NULL for anonymous.
#' @param capability Required capability.
#' @return One authorized dataset row; errors without revealing its existence.
#' @export
ribocrypt_access_resolve <- function(con, dataset, identity = NULL, capability = c("read", "download")) {
  access_scalar(dataset, "dataset")
  catalog <- ribocrypt_access_catalog(con, identity, match.arg(capability))
  result <- catalog[catalog$id == dataset, , drop = FALSE]
  if (nrow(result) != 1L) stop("Dataset unavailable or access denied.", call. = FALSE)
  result
}
