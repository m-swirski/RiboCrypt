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
                    paste("add_uorfs", input$add_uorfs, sep = "="),
                    paste("summary_track", input$summary_track, sep = "="),
                    paste("go", "TRUE", sep = "="),
                    sep = "&")

  page <- getPageFromURL(session)
  page <- ifelse(page == "", "", paste0("#", page))
  # Now combine
  url <- paste0(host, settings, page)
  print(paste("Copied text:", url))
  return(url)
}

clipboard_url_button <- function(input, session) {
  rclipButton(
    inputId = "clip",
    label = "Get URL",
    clipText = make_url_from_inputs(input, session),
    icon = icon("clipboard"),
    tooltip = "Get URL to share for this plot. Copied to clipboard (ctrl+v to paste)",
    placement = "top",
    options = list(delay = list(show = 800, hide = 100), trigger = "hover")
  )
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
  })
}

check_url_for_basic_parameters <- function() {
  with(rlang::caller_env(), {
    observeEvent(session$clientData$url_hash, {
      # Update experiment from url api
      page <- getPageFromURL(session)
      req(id == page || (page == "" && id == "browser") || (page == "MetaBrowser" && id == "browser_allsamp"))
      query <- getQueryString()
      tag <- "dff"
      value <- query[tag][[1]]

      if (is.null(input[[tag]]) || !is.null(value) && value != input[[tag]]
          && rv$exp != value) {
        print("Update experiment from url API")
        rv$exp <- value
      }
    }, priority = -5)

    observeEvent(session$clientData$url_hash, {
      page <- getPageFromURL(session)
      req(id == page || (page == "" && id == "browser") || (page == "MetaBrowser" && id == "browser_allsamp"))
      query <- getQueryString()
      print(paste("Page:", id))
      tag <- "gene"
      value <- query[tag][[1]]
      if (is.null(input[[tag]]) || !is.null(value) && value != input[[tag]]) {
        print(paste("Gene before:",isolate(input$gene)))
        print(paste("Update to:", value))
        gene_update_select(gene_name_list, selected = value)
        print(paste("Gene after:", isolate(input$gene)))
      }

      tag <- "tx"
      value <- query[tag][[1]]
      if (is.null(input[[tag]]) || !is.null(value) && value != input[[tag]]){
        # freezeReactiveValue(input, tag)
        tx_update_select(gene_name_list = gene_name_list, selected = value)
        print(isolate(input$gene))
      }
      tag <- "library"
      value <- query[tag][[1]]
      if (!is.null(value)){
        print(paste("Library update to:", value))
        value <- strsplit(x = value, ",")[[1]]
        if (length(value) > 0) {
          is_run_ids <- grep("SRR|ERR|DRR", value)
          l <- isolate(libs())
          matches_run <- matches_run_other <- TRUE
          if (length(is_run_ids) >  0) {
            print("Convert to ")
            run_ids <- runIDs(isolate(df()))
            matches_run <- run_ids %in% value[is_run_ids]
            matches_run_other <- value %in% run_ids
          }

          if (length(value) == 1 && value == "all") {
            value <- l
          } else {
            matches <- (l %in% value) | matches_run
            matches_other <- (value %in% l) | matches_run_other
            if (!all(matches)) {
              warning("Given libraries from URL are not part of this experiment:", paste(value[!matches_other], collapse = ", "))
              if (all(!matches_other)) {
                value <- l[1]
              } else value <- l[matches]
            }
          }
        }

        library_update_select(libs, selected = value)
        print(isolate(input$library))
      }

      tag <- "frames_type"
      value <- query[tag][[1]]
      if (!is.null(value)) {
        # freezeReactiveValue(input, tag)
        frame_type_update_select(value)
      }
      tag <- "kmer"
      value <- query[tag][[1]]
      if (!is.null(value)) {
        # freezeReactiveValue(input, tag)
        kmer_update_select(value)
      }
      tag <- "extendLeaders"
      value <- query[tag][[1]]
      if (!is.null(value)) {
        updateNumericInput(inputId = tag, value = value)
      }
      tag <- "extendTrailers"
      value <- query[tag][[1]]
      if (!is.null(value)) {
        updateNumericInput(inputId = tag, value = value)
      }
      tag <- "viewMode"
      value <- query[tag][[1]]
      if (!is.null(value)) {
        updateCheckboxInput(inputId = tag, value = as.logical(value))
      }
      tag <- "other_tx"
      value <- query[tag][[1]]
      if (!is.null(value)) {
        updateCheckboxInput(inputId = tag, value = as.logical(value))
      }
      tag <- "add_uorfs"
      value <- query[tag][[1]]
      if (!is.null(value)) {
        updateCheckboxInput(inputId = tag, value = as.logical(value))
      }
      tag <- "summary_track"
      value <- query[tag][[1]]
      if (!is.null(value)) {
        updateCheckboxInput(inputId = tag, value = as.logical(value))
      }

    }, priority = -10)

  })
}

check_url_for_go_on_init <- function() {
  with(rlang::caller_env(), {
    no_go_yet <- reactiveVal(TRUE)
    observeEvent(session$clientData$url_hash, {
      page <- getPageFromURL(session)
      req(id == page || (page == "" && id == "browser") || (page == "MetaBrowser" && id == "browser_allsamp"))
      query <- getQueryString()
      tag <- "go"
      value <- query[tag][[1]]
      if (!is.null(value)) {
        if (value[1] == TRUE) {
          print("Ready, set...")
          no_go_yet(FALSE)
          browser_options["plot_on_start"] <- "FALSE"
          print("Set plot_on_start to FALSE")
        }
      }
    }, ignoreNULL = TRUE, ignoreInit = FALSE, priority = -100)
    # Timer for running plot, we have to wait for setup to finish
    rtimer <- reactiveTimer(1000)
    timer <- reactive({req(no_go_yet() == FALSE);print("Timer activated!"); rtimer()}) %>% bindEvent(rtimer(), ignoreInit = TRUE)

    observeEvent(timer(), {
      if (!no_go_yet()) {
        req(input$gene != "")
        print(paste("Fire gene: ", isolate(input$gene)))
        query <- getQueryString()
        tag <- "gene"
        value <- query[tag][[1]]
        if (!is.null(value)) req(input$gene == value)
        req(input$tx != "" && !is.null(input$tx))
        tag <- "tx"
        value <- query[tag][[1]]
        if (!is.null(value)) req(input$tx == value)
        print(paste("Fire tx: ", isolate(input$tx)))
        print("Fire button!")
        shinyjs::click("go")
        no_go_yet(TRUE)
      }
    }, ignoreInit = TRUE, ignoreNULL = TRUE, priority = -200)

  })
}
