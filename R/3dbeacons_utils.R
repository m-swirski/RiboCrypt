beacon_main_api_url <- function() {
  'https://www.ebi.ac.uk/pdbe/pdbe-kb/3dbeacons/api/uniprot/summary/'
}

#' Fetch summary of uniprot id
#' @param qualifier uniprot ids
#' @param provider "pdbe", alternatives: "alphafold", "all"
#' @return a character of json
fetch_summary <- function(qualifier, provider = "pdbe") {
  provider_url <-
  if (provider == "all") {
    stop("This does not work yet!")
    ""
  } else paste0("?provider=", provider)

  stringr::str_interp(
    paste0(beacon_main_api_url(), "${qualifier}", ".json", provider_url),
    list(qualifier = qualifier)
  ) %>% RCurl::getURI() %>% jsonlite::fromJSON()
}

model_identifiers_from_summary <- function(data) {
  data$model_identifiers
}

model_categories_from_summary <- function(data) {
  data$model_category
}

model_urls_from_summary <- function(data) {
  data$structures$summary$model_url
}

model_formats_from_summary <- function(data) {
  data$model_format
}

entity_identifiers_from_summary <- function(data) {
  mapply(function (x) {
    x$identifier
  },
  data$structures$summary$entities)
}

entity_identifier_categories_from_summary <- function(data) {
  mapply(function (x) {
    x$identifier_category
  },
  data$structures$summary$entities)
}

chain_ids_from_summary <- function(data) {
  mapply(function (x) {
    x$chain_ids
  },
  data$structures$summary$entities)
}
