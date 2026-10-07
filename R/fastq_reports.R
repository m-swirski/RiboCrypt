rc_fastq_report_response <- function(path, directory, context) {
  tryCatch({
    rc_access_check(context)
    root <- normalizePath(directory, mustWork = TRUE)
    file <- normalizePath(path, mustWork = TRUE)
    if (!startsWith(file, paste0(root, "/")) || dir.exists(file))
      stop("Report unavailable or access denied.", call. = FALSE)
    list(status = 200L, headers = list(
      "Content-Type" = "text/html; charset=UTF-8", "Cache-Control" = "no-store",
      "Content-Security-Policy" = "sandbox allow-scripts", "X-Content-Type-Options" = "nosniff"),
      body = readBin(file, "raw", n = file.info(file)$size))
  }, error = function(e) rc_access_response())
}

rc_fastq_report_url <- function(path, directory, session = shiny::getDefaultReactiveDomain()) {
  if (is.null(session)) stop("A report requires an active session.", call. = FALSE)
  context <- session$userData$ribocrypt_access
  # Register only this report; never expose its directory through a static resource.
  session$registerDataObj("fastq_report", NULL, function(data, req) {
    rc_fastq_report_response(path, directory, context)
  })
}
