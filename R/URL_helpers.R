getPageFromURL <- function(session = NULL, url = session$clientData$url_hash) {
  utils::URLdecode(sub("#", "", session$clientData$url_hash))
}
