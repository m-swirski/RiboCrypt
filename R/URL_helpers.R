getPageFromURL <- function(session = NULL, url = session$clientData$url_hash) {
  utils::URLdecode(sub("#", "", session$clientData$url_hash))
}

make_url_from_inputs <- function(input, session) {
  host <- session$clientData$url_hostname
  if (host != "ribocrypt.org") {
    host <- paste0(host, ":", session$clientData$url_port)
  }
  settings <- "/?"
  settings <- paste(settings,
                    paste("dff", input$dff, sep = "="),
                    paste("gene", input$gene, sep = "="),
                    paste("tx", input$tx, sep = "="),
                    paste("library", paste(input$library, collapse = ","), sep = "="),
                    paste("frames_type", input$frames_type, sep = "="),
                    paste("kmer", input$kmer, sep = "="),
                    paste("extendLeaders", input$extendLeaders, sep = "="),
                    paste("extendTrailers", input$kmer, sep = "="),
                    paste("viewMode", input$viewMode, sep = "="),
                    paste("other_tx", input$other_tx, sep = "="),
                    paste("go", "TRUE", sep = "="),
                    sep = "&")

  page <- getPageFromURL(session)
  page <- ifelse(page == "", "", paste0("#", page))
  # Now combine
  url <- paste0(host, settings, page)
  print(paste("Copied text:", url))
  return(url)
}

#' Make the URL field reactive to page given
#'
#' Currently does not support update of input fields other than selected page
#' @noRd
reactive_url <- function() {
  with(rlang::caller_env(), {
    observeEvent(session$clientData$url_hash, {
      currentHash <- getPageFromURL(session)
      if (is.null(input$navbarID) || !is.null(currentHash) && currentHash != input$navbarID){
        freezeReactiveValue(input, "navbarID")
        updateNavbarPage(session, "navbarID", selected = currentHash)
      }
    }, priority = 1)
    observeEvent(input$navbarID, {
      currentHash <- getPageFromURL(session)
      pushQueryString <- paste0("#", input$navbarID)
      if(is.null(currentHash) || currentHash != input$navbarID){
        freezeReactiveValue(input, "navbarID")
        updateQueryString("?", mode = "replace", session)
        updateQueryString(pushQueryString, mode = "push", session)
      }
    }, priority = 0, ignoreInit = TRUE)
  }
  )
}
