access_secret_equal <- function(actual, expected) {
  if (!is.character(actual) || length(actual) != 1L || is.na(actual)) return(FALSE)
  a <- charToRaw(actual)
  b <- charToRaw(expected)
  if (length(a) != length(b)) return(FALSE)
  difference <- 0L
  for (i in seq_along(b)) difference <- bitwOr(difference, bitwXor(as.integer(a[i]), as.integer(b[i])))
  difference == 0L
}

#' Verify Identity Supplied by a Trusted Authentication Gateway
#' @param request Shiny/httpuv request containing HTTP_RIBOCRYPT_GATEWAY_SECRET
#'   and optionally HTTP_RIBOCRYPT_SUBJECT. NGINX must overwrite both headers.
#' @param issuer Fixed OIDC issuer configured by the administrator.
#' @param gateway_secret Deployment secret of at least 32 bytes. Never expose it
#'   in browser code, URLs or logs. Restrict backend access to the proxy.
#' @return NULL for anonymous requests, otherwise an issuer/subject list.
#' @details This verifies the gateway, not an OIDC token. The gateway must
#'   authenticate the user and supply the immutable OIDC subject. Do not trust
#'   email, URL parameters or headers arriving directly from a browser.
#' @export
ribocrypt_access_identity <- function(request, issuer, gateway_secret) {
  access_scalar(issuer, "issuer")
  access_scalar(gateway_secret, "gateway_secret")
  if (nchar(gateway_secret, type = "bytes") < 32L)
    stop("Gateway secret must contain at least 32 bytes.", call. = FALSE)
  if (!access_secret_equal(request$HTTP_RIBOCRYPT_GATEWAY_SECRET, gateway_secret))
    stop("Untrusted authentication gateway.", call. = FALSE)
  subject <- request$HTTP_RIBOCRYPT_SUBJECT
  if (is.null(subject) || identical(subject, "")) return(NULL)
  list(issuer = issuer, subject = access_scalar(subject, "subject"))
}
