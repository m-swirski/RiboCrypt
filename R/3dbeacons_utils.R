beacon_main_api_url <- function() {
  'https://www.ebi.ac.uk/pdbe/pdbe-kb/3dbeacons/api/uniprot/summary/'
}

alphafold_main_dir <- function() "https://alphafold.ebi.ac.uk/files/"

get_protein_structure <- function(gene_name_list, selectedRegion, selectedTX, region_dir,
                                  protein_structure_dir, df, customRegions) {
  print("Selecting local or online pdb")
  uniprot_id <- gene_name_list[gene_name_list$value == selectedRegion]$uniprot_id

  # TODO wrap this in a general fetches
  protein_clicked_info(selectedRegion, selectedTX, uniprot_id)
  online <- isTruthy(uniprot_id)
  pdb_files <- NULL
  if (online) {
    pdb_files <- try(get_protein_structure_online(uniprot_id), silent = TRUE)
  }

  online <- isTruthy(pdb_files)
  if (online) {
    res <- pdb_files
    print("Structure searched online")
  } else {
    res <- get_protein_structure_local_on_disk(selectedRegion,
                                               region_dir, protein_structure_dir,
                                               df, customRegions)
    print("Structure searched local")
  }
  print(res)

  # if (is.null(input$structureViewerSelector)) {
  #   head(structure_variants())
  # } else input$structureViewerSelector
  print(paste("Selected variant: ", head(res, 1)))
  attr(res, "selected") <- head(res, 1)
  return(res)
}

#' Download pdb by uniprot id
get_protein_structure_online <- function(uniprot_id) {
  results <- get_protein_structure_urls(uniprot_id)
  tempfiles <- make_tempfile_paths_for_protein_structures(results)
  download_protein_structures_by_url(tempfiles, results)
  return(tempfiles)
}

get_protein_structure_urls <- function(uniprot_id) {
  print("Fetching protein structure URLs..")
  model_urls <-
    fetch_summary(uniprot_id) %>%
    model_urls_from_summary()
  # Convert alphafold cif to pdb (NGLviewR does not work with cif there)
  alphafold_id <- grep(alphafold_main_dir(), model_urls)
  if (length(alphafold_id) > 0) {
    model_urls[alphafold_id] <- gsub("cif$", "pdb", model_urls[alphafold_id])
  }
  return(model_urls)
}

#' Fetch summary of uniprot id from provider
#' @param qualifier uniprot ids
#' @param provider "pdbe", alternatives: "alphafold", "all"
#' @return a character of json
#' @importFrom RCurl getURI
#' @importFrom jsonlite fromJSON
fetch_summary <- function(qualifier, provider = "alphafold") {
  provider_url <-
  if (provider == "all") {
    stop("This does not work yet!")
    ""
  } else paste0("?provider=", provider)
  summary <- stringr::str_interp(
    paste0(beacon_main_api_url(), "${qualifier}", ".json", provider_url),
    list(qualifier = qualifier)
  ) %>% RCurl::getURI()
  error_codes <- "502 Bad gateway|503 Service Temporarily Unavailable"
  if (grepl(error_codes, summary, ignore.case = TRUE)) {
    warning("Protein structure Provider webpage is down for the moment!")
    summary <- "{}"
  }

  jsonlite::fromJSON(summary)
}

make_tempfile_paths_for_protein_structures <- function(results) {
  model_labels <- mapply(
    function(x) {
      str_sub(x, start = gregexpr("/", x) %>% unlist() %>% last() + 1)
    },
    results
  )
  model_paths <- rep(tempfile(pattern = "structure"), length(model_labels))
  model_paths <- paste0(model_paths, ".", file_ext(results))
  result <- model_paths
  names(result) <- model_labels
  return(result)
}

download_protein_structures_by_url <- function(tempfiles, results) {
  print("Fetching structures..")
  tmp_paths <- unname(tempfiles)
  names(tmp_paths) <- results
  mapply(
    function(x) {
      httr::GET(x, httr::write_disk(tmp_paths[x], overwrite = TRUE))
    },
    results
  )
}

get_protein_structure_local_on_disk <- function(selectedRegion, region_dir,
                                                protein_structure_dir, df,
                                                customRegions) {
  paths <- character()
  if (!isTruthy(selectedRegion)) return(paths)
  pdb_input <- grepl("\\.pdb$", selectedRegion)
  uorf_clicked <- length(grep("^U[0-9]+$", selectedRegion)) == 1
  translon_clicked <- length(grep("^T[0-9]+$", selectedRegion)) == 1
  if (pdb_input) {
    paths <- selectedRegion
  } else if (uorf_clicked) {
    paths <- list.files(region_dir, full.names = TRUE)
    paths <- paths[grep(paste0("^uorf_",selectedRegion, ".pdb$"), basename(paths))]
    if (length(paths) == 0) warning("No local protein structure for this uORF!")
  } else if (translon_clicked) {
    linker_file <- file.path(refFolder(df), "predicted_translons",
                             "predicted_translons_with_sequence_pep_linker.fst")
    if (!file.exists(linker_file)) {
      print("No translon peptide linker file for organism!")
      return("")
    }
    selected_grl <- customRegions
    selected_as_coord <- as.character(selected_grl[names(selected_grl) == selectedRegion])
    linker_dt <- fst::read_fst(linker_file, as.data.table = TRUE)

    paths <- coordinates_to_pep_id_path(selected_as_coord, linker_dt, protein_structure_dir)
  } else {
    paths <- paths[-grep("uorf", paths)] # Remove uorf structures
  }

  path_labels <- mapply(
    function(x) {
      str_sub(x, start = gregexpr("/", x) %>% unlist() %>% last() + 1)
    },
    basename(paths)
  )
  result <- paths
  names(result) <- path_labels
  result
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
