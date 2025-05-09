#' Make URL from shiny reactive input
#' @noRd
make_url_from_inputs <- function(input, session) {
  host <- getHostFromURL(session)
  parameters <- make_url_from_inputs_parameters(input)
  page <- getPageFromURL(session, with_hash = TRUE)

  # Now combine
  url <- paste0(host, parameters, page)
  print(paste("URL:", url))
  return(url)
}

getPageFromURL <- function(session = NULL, url = session$clientData$url_hash,
                           with_hash = FALSE) {
  page <- utils::URLdecode(sub("#", "", session$clientData$url_hash))
  if (with_hash) page <- ifelse(page == "", "", paste0("#", page))
  return(page)
}

getHostFromURL <- function(session) {
  host <- session$clientData$url_hostname
  if (host == "ribocrypt.org") {
    host <- paste0("https://", host)
  } else if (host == "ribocrypt.neutrino.re") {
    host <- paste0("https://", host, "/app/ribocrypt")
  } else { # Else local user / other server
    host <- paste0("http://", host, ":", session$clientData$url_port)
  }
  return(host)
}

make_url_from_inputs_parameters <-function(input, go = TRUE, settings = "/?") {
  paste0(settings,
  paste(paste("dff", input$dff, sep = "="),
        paste("gene", input$gene, sep = "="),
        paste("tx", input$tx, sep = "="),
        paste("library", paste(input$library, collapse = ","), sep = "="),
        paste("frames_type", input$frames_type, sep = "="),
        paste("kmer", input$kmer, sep = "="),
        paste("log_scale", input$log_scale, sep = "="),
        paste("log_scale_protein", input$log_scale_protein, sep = "="),
        paste("extendLeaders", input$extendLeaders, sep = "="),
        paste("extendTrailers", input$extendTrailers, sep = "="),
        paste("viewMode", input$viewMode, sep = "="),
        paste("other_tx", input$other_tx, sep = "="),
        paste("add_uorfs", input$add_uorfs, sep = "="),
        paste("add_translon", input$add_translon, sep = "="),
        paste("genomic_region", sub("\\+;", "p;", sub("\\+$", "p", input$genomic_region)), sep = "="),
        paste("zoom_range", sub("\\+$", "p", input$zoom_range), sep = "="),
        paste("customSequence", input$customSequence, sep = "="),
        paste("phyloP", input$phyloP, sep = "="),
        paste("collapsed_introns_width", input$collapsed_introns_width, sep = "="),
        paste("collapsed_introns", input$collapsed_introns, sep = "="),
        paste("go", go, sep = "="),
        sep = "&"))
}

clipboard_url_button <- function(input, session) {
  rclipButton(
    inputId = "clip",
    label = "Get URL",
    clipText = make_url_from_inputs(input, session),
    icon = icon("clipboard"),
    tooltip = "Get URL to share for this plot. Copied to clipboard (ctrl+v to paste)",
    placement = "top",
    options = list(delay = list(show = 600, hide = 100), trigger = "hover")
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

    check_go_flag()
  })
}

check_url_for_basic_parameters <- function() {
  with(rlang::caller_env(), {

    observeEvent(session$clientData$url_hash, { # PAGE: Browsers
      # Update experiment from url api
      page <- getPageFromURL(session)
      req(id == page || (page == "" && id == "browser") || (page == "MegaBrowser" && id == "browser_allsamp"))
      query <- getQueryString()

      tag <- "dff"
      value <- query[tag][[1]]
      if (is.null(input[[tag]]) || !is.null(value) && value != input[[tag]]
          && rv$exp != value) {
        print("Update experiment from url API")
        rv$exp <- value
      }
    }, priority = -5)

    observeEvent(session$clientData$url_hash, { # PAGE: Predicted translons tables
      # Update experiment from url api
      page <- getPageFromURL(session)
      req(id == "predicted_translons" && page %in% c("predicted_translons", "Predicted Translons"))
      query <- getQueryString()
      req(length(query) > 0)
      message("URL with predicted translons!")
      tag <- "dff"
      value <- query[tag][[1]]
      if (is.null(input[[tag]]) || !is.null(value) && value != input[[tag]]) {
        updateSelectizeInput(
          inputId = "dff",
          selected = value,
          server = FALSE
        )
      }
    }, priority = -6)

    observeEvent(session$clientData$url_hash, { # PAGE: Predicted translons tables
      # Update experiment from url api
      page <- getPageFromURL(session)
      req(id == "predicted_translons" && page %in% c("predicted_translons", "Predicted Translons"))
      query <- getQueryString()
      req(length(query) > 0)

      tag <- "download_format"
      value <- query[tag][[1]]
      if (is.null(input[[tag]]) || !is.null(value) && value %in% c("xlsx", "csv")) {
        type <- ifelse(value == "xlsx", "excel", value)
        trigger_input <- paste0("trigger_download_", type)
        download_button <- paste0("download_", type)
        message("Downloading translon dataset!")
        download_trigger(type)
      }
    }, priority = -150)

    observeEvent(session$clientData$url_hash, { # PAGE: Browsers
      page <- getPageFromURL(session)
      req(id == page || (page == "" && id == "browser") || (page == "MegaBrowser" && id == "browser_allsamp"))
      query <- getQueryString()
      req(length(query) > 0)
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
      if (!is.null(value)) {
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
        frame_type_update_select(value)
      }
      tag <- "kmer"
      value <- query[tag][[1]]
      if (!is.null(value)) {
        kmer_update_select(value)
      }

      # Numeric box updates
      for (tag in c("extendLeaders", "extendTrailers", "collapsed_introns_width")) {
        value <- query[tag][[1]]
        if (!is.null(value)) {
          updateNumericInput(inputId = tag, value = value)
        }
      }
      # Free character box updated
      for (tag in c("customSequence", "genomic_region", "zoom_range")) {
        value <- query[tag][[1]]
        if (!is.null(value)) {
          if (tag %in% c("genomic_region", "zoom_range")) value <-  sub("p;", "+;", sub("p$", "+", value))
          updateTextInput(inputId = tag, value = as.character(value))
        }
      }

      # Checkbox updates
      for (tag in c("viewMode", "other_tx", "add_uorfs", "add_translon","summary_track",
                    "log_scale", "log_scale_protein","phyloP", "collapsed_introns")) {
        value <- query[tag][[1]]
        if (!is.null(value)) {
          updateCheckboxInput(inputId = tag, value = as.logical(value))
        }
      }
    }, priority = -10)

  })
}

check_url_for_go_on_init <- function() {
  with(rlang::caller_env(), {
    no_go_yet <- reactiveVal(TRUE)
    observeEvent(session$clientData$url_hash, {
      page <- getPageFromURL(session)
      req(id == page || (page == "" && id == "browser") || (page == "MegaBrowser" && id == "browser_allsamp"))
      query <- getQueryString()
      tag <- "go"
      value <- query[tag][[1]]
      if (!is.null(value)) {
        if (value[1] == TRUE) {
          print("Ready, set...")
          no_go_yet(FALSE)
          browser_options["plot_on_start"] <- "FALSE"
        }
      }
    }, ignoreNULL = TRUE, ignoreInit = FALSE, priority = -100)
    # Timer for running plot, we have to wait for setup to finish
    rtimer <- reactiveTimer(1000)
    timer <- reactive({req(no_go_yet() == FALSE);print("Timer activated!"); rtimer()}) %>%
      bindEvent(rtimer(), ignoreInit = TRUE)

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

check_go_flag <- function() {
  with(rlang::caller_env(), {
    go_flag <- isolate(reactive({
      query <- getQueryString()
      tag <- "go"
      as.logical(query[tag][[1]][1])
    })())

    if (length(go_flag) > 0 && go_flag %in% c(TRUE, FALSE)) {
      over_ride_plot_on_start <- as.logical(go_flag) == FALSE && as.logical(browser_options["plot_on_start"])
      if (over_ride_plot_on_start) {
        browser_options["plot_on_start"] <- "FALSE"
        print("Set plot_on_start to FALSE")
      }
    }

  })
}

#' Browse a gene on Ribocrypt webpage
#'
#' Can also disply local RiboCrypt app
#' @inheritParams make_rc_url
#' @param browser getOption("browser")
#' @return browseURL, opens browse with page
#' @export
#' @examples
#' browseRC("ATF4", "ENSG00000128272")
#'
browseRC <- function(symbol = NULL, gene_id = NULL, tx_id = NULL,
                     exp = "all_merged-Homo_sapiens_modalities",
                     libraries = NULL, leader_extension = 0, trailer_extension = 0,
                     viewMode = FALSE, other_tx = FALSE,
                     plot_on_start = TRUE, frames_type = "columns", kmer=1,
                     host = "https://ribocrypt.org",
                     browser = getOption("browser")) {

  full_url <- make_rc_url(symbol, gene_id, tx_id, exp, libraries,
                          leader_extension, trailer_extension,
                          viewMode, other_tx, plot_on_start, frames_type,
                          kmer, host)
  browseURL(full_url, browser = browser)
}


#' Create URL to browse a gene on Ribocrypt webpage
#'
#' Can also make url for local RiboCrypt app'
#' On the actuall app, the function make_url_from_inputs is used on
#' the shiny reactive input object. This one is for manual use.
#' @inheritParams multiOmicsPlot_list
#' @param symbol gene symbol, default NULL
#' @param gene_id gene symbol, default NULL
#' @param tx_id gene symbol, default NULL
#' @param exp experiment name, default "all_merged-Homo_sapiens_modalities"
#' @param libraries NULL, default to first in experiment, c("RFP","RNA") would add RNA to default.
#' @param plot_on_start logical, default TRUE. Plot gene when opening browser.
#' @param frames_type "columns"
#' @param viewMode FALSE (transcript view), TRUE gives genomic.
#' @param other_tx FALSE, show all other annotation in region (isoforms etc.)
#' @param kmer integer, default 1 (no binning), binning size of windows, to smear out the signal.
#' @param add_translons logical, default FALSE. If TRUE, add translons predicted sequences in annotation.
#' @param zoom_range character, zoom values.
#' @param host url, default "https://ribocrypt.org". Set to localhost for local version.
#' @return character, URL.
#' @export
#' @examples
#' make_rc_url("ATF4", "ENSG00000128272")
make_rc_url <- function(symbol = NULL, gene_id = NULL, tx_id = NULL,
                        exp = "all_merged-Homo_sapiens_modalities",
                        libraries = NULL, leader_extension = 0, trailer_extension = 0,
                        viewMode = FALSE, other_tx = FALSE,
                        plot_on_start = TRUE, frames_type = "columns", kmer=1,
                        add_translons = FALSE, zoom_range = NULL,
                        host = "https://ribocrypt.org") {
  # TODO: Update to use input function above, to avoid writing twice
  if (any(is.null(symbol) & is.null(gene_id)))
    stop("At least on of symbol and gene_id must be defined!")
  settings <- "/?"
  settings <- paste0(settings,
                    paste(paste("dff", exp, sep = "="),
                    paste("frames_type", frames_type, sep = "="),
                    paste("kmer", kmer, sep = "="),
                    paste("extendLeaders", leader_extension, sep = "="),
                    paste("extendTrailers", trailer_extension, sep = "="),
                    paste("viewMode", viewMode, sep = "="),
                    paste("other_tx", other_tx, sep = "="),
                    paste("add_translon", add_translons, sep = "="),
                    sep = "&"))
  if (!is.null(zoom_range)){
    zoom_range <- sub("\\+;", "p;", sub("\\+$", "p", zoom_range))
    settings <- paste(settings, paste("zoom_range", zoom_range, sep = "="), sep = "&")
  }

  prefix_url <- paste0(host, settings)
  if (!is.null(symbol)) {
    dt <- data.table(gene_id, symbol)
    dt[, merged_name := do.call(paste, .SD, ), .SDcols = c(2,1)]
    dt[, merged_name := sub(" ",  "-", merged_name, fixed = TRUE)]
    dt[, merged_name := sub("(^-)|(^NA-)",  "", merged_name, perl = TRUE)]
    gene <- dt$merged_name
  } else {
    gene <- gene_id
  }

  if (all(!is.null(tx_id))) tx_id <- paste0("&tx=", tx_id)
  if (all(!is.null(libraries))) libraries <- paste0("&library=", paste(libraries, collapse = ","))

  select <- paste0("&gene=", gene, tx_id, libraries)
  plot_on_start <- paste0("&go=", as.logical(plot_on_start))
  full_url <- paste0(prefix_url, select, plot_on_start)
  return(full_url)
}
