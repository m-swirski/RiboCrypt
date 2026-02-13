#' Make URL from shiny reactive input
#' @noRd
make_url_from_inputs <- function(input, session) {
  host <- getHostFromURL(session)
  parameters <- make_url_from_inputs_parameters(input)
  page <- getPageFromURL(session, with_hash = TRUE)

  # Now combine
  url <- paste0(host, parameters, page)
  cat(paste("URL:", url), sep = "\n")
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
    port <- session$clientData$url_port
    pathname <- sub("/$", "", session$clientData$url_pathname)
    host <- paste0("http://", host, ":", port, pathname)
  }
  return(host)
}

make_url_from_inputs_parameters <-function(input, go = TRUE, settings = "/?") {
  paste0(settings,
  paste(paste("dff", input$dff, sep = "="),
        paste("gene", input$gene, sep = "="),
        paste("tx", input$tx, sep = "="),
        paste("library", paste(input$library, collapse = ","), sep = "="),
        paste("unique_align", paste(input$unique_align, collapse = ","), sep = "="),
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
        paste("add_translons_transcode", input$add_translons_transcode, sep = "="),
        paste("genomic_region", sub("\\+;", "p;", sub("\\+$", "p", input$genomic_region)), sep = "="),
        paste("zoom_range", sub("\\+$", "p", input$zoom_range), sep = "="),
        paste("customSequence", input$customSequence, sep = "="),
        paste("phyloP", input$phyloP, sep = "="),
        paste("mapability", input$mapability, sep = "="),
        paste("colors", input$colors, sep = "="),
        paste("summary_track", input$summary_track, sep = "="),
        paste("summary_track_type", input$summary_track_type, sep = "="),
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
    cat("Server startup: "); print(round(Sys.time() - time_before, 2))
    this_time_before <- Sys.time()

    # Page routing API
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
    # Argument API
    check_url_for_basic_parameters()
  })
}

check_url_for_basic_parameters <- function() {
  with(rlang::caller_env(), {

  url_args <- isolate(getQueryString())
  if (length(url_args) > 0) {
    query <- url_args

    tag <- "dff"
    value <- query[tag][[1]]
    if (!is.null(value)) {
      print(paste("Update experiment from url API"))
      if (!(value %in% all_exp$name)) stop("You specified invalid experiment from URL API!")
      browser_options["default_experiment"] <- value
    }
    message("From URL:", browser_options["default_experiment"])

    tag <- "gene"
    value <- query[tag][[1]]
    if (!is.null(value)) {
      print(paste("Gene before:", browser_options["default_gene"]))
      print(paste("Update to:", value))
      browser_options["default_gene"] <- value
    }

    tag <- "tx"
    value <- query[tag][[1]]
    if (!is.null(value)){
      print(paste("Tx before:", browser_options["default_isoform"]))
      print(paste("Update to:", value))
      browser_options["default_isoform"] <- value
    }
    tag <- "library"
    value <- query[tag][[1]]
    if (!is.null(value)) {
      print(paste("Library update to:", value))
      browser_options["default_libs"] <- value
    }

    tag <- "go"
    value <- query[tag][[1]]
    if (!is.null(value)) {
      temp_value <- as.logical(value)
      if (!is.na(temp_value)) {
        browser_options["plot_on_start"] <- as.character(value)
      }
    }

    tag <- "search"
    value <- query[tag][[1]]
    if (!is.null(value)) {
      temp_value <- as.character(value)
      if (!is.na(temp_value)) {
        browser_options["search_on_init"] <- as.character(value)
      }
    }

  }
  })
}

browser_specific_url_checker <- function() {
  with(rlang::caller_env(), {
    # url tags
    if (id == "browser") {
      url_args <- isolate(getQueryString())
      if (length(url_args) > 0) {
        query <- url_args
        print("Browser specific URL check")
        for (tag in c("frames_type", "summary_track_type")) {
          value <- query[tag][[1]]
          if (!is.null(value)) {
            frame_type_update_select(value, tag)
          }
        }
        for (tag in c("colors")) {
          value <- query[tag][[1]]
          if (!is.null(value)) {
            updateSelectizeInput(inputId = tag, selected = value)
          }
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
            if (tag %in% c("genomic_region", "zoom_range")) value <-  gsub("p;", "+;", sub("p$", "+", value))
            updateTextInput(inputId = tag, value = as.character(value))
          }
        }

        # Checkbox updates
        for (tag in c("viewMode", "other_tx", "add_uorfs", "add_translon",
                      "add_translons_transcode", "summary_track",
                      "log_scale", "log_scale_protein","phyloP", "mapability",
                      "collapsed_introns", "unique_align")) {
          value <- query[tag][[1]]
          if (!is.null(value)) {
            updateCheckboxInput(inputId = tag, value = as.logical(value))
          }
        }
      }
    }

    # User info is checked through browser
    if (id == "browser") {
      kickoff <- reactiveVal(FALSE)
      fired <- reactiveVal(FALSE)
      observeEvent(list(input$gene, input$tx, input$library),
                   go_when_input_is_ready(input, browser_options, fired, kickoff, libs),
                   ignoreInit = TRUE, ignoreNULL = TRUE)

      user_info <- reactive({
        is_cellphone <- grepl("Android|iPhone|iPad|iPod",
                              input$js_user_agent,
                              ignore.case = TRUE)
        list(
          userAgent = input$js_user_agent,
          is_cellphone = is_cellphone,
          width     = session$clientData$`output_browser-c_width`,
          height    = session$clientData$`output_browser-c_height`
        )
      })

      # Optional: print when something changes
      observe({
        cat("User info updated:\n")
        str(user_info())
      })
    }
  })
}

translon_specific_url_checker <- function() {
  with(rlang::caller_env(), {
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
        download_trigger(type) # Defined outside
      }
      tag <- "dff"
      value <- query[tag][[1]]
      if (!is.null(value)) {
        print(paste("Update experiment from url API"))
        if (!(value %in% all_exp$name)) stop("You specified invalid experiment from URL API!")
        browser_options["default_experiment_translon"] <- value
      }

    }, priority = -150)
  })
}



libraries_string_split <- function(value, libs) {
  value <- strsplit(x = value, ",")[[1]]
  if (length(value) > 0 && !anyNA(value)) {
    is_run_ids <- grep("SRR|ERR|DRR", value)
    l <- libs
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
  } else value <- libs[1]
  if (!all(value %in% libs))
    stop("You defined libraries to use, but some of those are not valid names,",
         " in the selected experiment!")
  return(value)
}

#' Browse a gene on Ribocrypt webpage
#'
#' Can also disply local RiboCrypt app if specified in the 'host' argument.
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
                     add_translons = FALSE, add_translons_transcode = FALSE,
                     zoom_range = NULL,
                     genomic_region = NULL,
                     host = "https://ribocrypt.org",
                     browser = getOption("browser")) {

  full_url <- make_rc_url(symbol, gene_id, tx_id, exp, libraries,
                          leader_extension, trailer_extension,
                          viewMode, other_tx, plot_on_start, frames_type,
                          kmer, add_translons, add_translons_transcode,
                          zoom_range, genomic_region, host)
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
#' @param genomic_region character, region to display
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
                        add_translons = FALSE, add_translons_transcode = FALSE,
                        zoom_range = NULL,
                        genomic_region = NULL,
                        host = "https://ribocrypt.org") {
  # TODO: Update to use input function above, to avoid writing twice
  if (any(is.null(symbol) & is.null(gene_id) & is.null(genomic_region)))
    stop("At least on of symbol, gene_id and genomic_region must be defined!")
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
                    paste("add_translons_transcode", add_translons, sep = "="),
                    sep = "&"))
  if (!is.null(zoom_range)){
    zoom_range <- gsub("\\+;", "p;", sub("\\+$", "p", zoom_range))
    settings <- paste(settings, paste("zoom_range", zoom_range, sep = "="), sep = "&")
  }
  if (!is.null(genomic_region)){
    genomic_region <- gsub("\\+;", "p;", sub("\\+$", "p", genomic_region))
    settings <- paste(settings, paste("genomic_region", genomic_region, sep = "="), sep = "&")
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
  if (all(!is.null(gene))) gene_id <- paste0("&gene=", gene)
  select <- paste0(gene_id, tx_id, libraries)
  plot_on_start <- paste0("&go=", as.logical(plot_on_start))
  full_url <- paste0(prefix_url, select, plot_on_start)
  return(full_url)
}
