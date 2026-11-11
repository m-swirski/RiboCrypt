# Optional landing page. The experiment and metadata catalogs always come from
# the current session, including access-control filtering.
rc_home_host <- function(x) sub(" x .*$", "", as.character(x))

rc_home_groups <- function(x) {
  out <- rep("Single experiments", length(x))
  out[grepl("^all_merged-", x)] <- "Merged tracks"
  out[grepl("^all_merged-.*_modalities$", x)] <- "Modalities"
  out[grepl("^all_samples-", x)] <- "Collections"
  out[grepl("^all_(cell_lines|tissues|inhibitors|translation_phase|sample_origin|cancer_status)", x)] <- "Collections"
  out
}

rc_home_search <- function(m, q) {
  if (is.null(m) || !nrow(m) || !"Run" %in% names(m) || nchar(trimws(q)) < 2L) return(character())
  q <- tolower(trimws(q)); cols <- names(m)[vapply(m,function(x) is.character(x)||is.factor(x),logical(1))]
  keep <- rep(FALSE,nrow(m))
  for (col in cols) keep <- keep | (!is.na(m[[col]]) & grepl(q,tolower(as.character(m[[col]])),fixed=TRUE))
  unique(m$Run[keep & !is.na(m$Run)])
}

rc_home_catalog <- function(all_exp, metadata) {
  sp <- unique(rc_home_host(all_exp$organism)); sp <- sp[!is.na(sp) & nzchar(sp)]
  taxonomy <- utils::read.csv(system.file("home", "species_clades.csv", package="RiboCrypt"), stringsAsFactors=FALSE)
  out <- data.frame(name=sp, clade=taxonomy$clade[match(sp,taxonomy$species)], libraries=0L, studies=0L, stringsAsFactors=FALSE)
  out$clade[is.na(out$clade)] <- "Unclassified"
  if (!is.null(metadata) && "ScientificName" %in% names(metadata)) {
    m <- data.table::copy(metadata); m[, .home_species := rc_home_host(ScientificName)]
    stats <- m[, .(libraries=if ("Run" %in% names(m)) data.table::uniqueN(Run, na.rm=TRUE) else .N,
      studies=if ("study_accession" %in% names(m)) data.table::uniqueN(study_accession, na.rm=TRUE) else NA_integer_), by=.home_species]
    i <- match(out$name,stats$.home_species); ok <- !is.na(i)
    out$libraries[ok] <- stats$libraries[i[ok]]; out$studies[ok] <- stats$studies[i[ok]]
  }
  out
}

rc_home_ui <- function(home_html, all_exp, metadata) {
  shiny::tabPanel("Home", value="home", icon=shiny::icon("house"),
    shiny::tags$head(shiny::tags$link(rel="stylesheet",href="rc-home/home.css"),
                     shiny::tags$script(src="rc-home/home.js")),
    shiny::div(class="rc-home", shiny::includeHTML(home_html),
      shiny::fluidRow(
        shiny::column(7, shiny::div(class="rc-home-card", shiny::h3("A world of translation"),
          shiny::p("Choose a species on the wheel or use the search."),
          shiny::tags$svg(id="rc-home-wheel", viewBox="0 0 620 620", role="group", `aria-label`="Species wheel"),
          shiny::div(id="rc-home-tooltip",class="rc-home-tooltip"), shiny::div(id="rc-home-legend",class="rc-home-legend"))),
        shiny::column(5, shiny::div(class="rc-home-card", shiny::h3("Choose your data"),
          shiny::selectizeInput("home_species", "Species", choices=sort(unique(rc_home_host(all_exp$organism))),
            selected=if ("Homo sapiens" %in% rc_home_host(all_exp$organism)) "Homo sapiens" else rc_home_host(all_exp$organism)[1],
            width="100%", options=list(placeholder="Search species")),
          shiny::uiOutput("home_species_info"),
          shiny::selectizeInput("home_experiment", "Experiment", choices=character(), width="100%", options=list(maxOptions=100)),
          shiny::uiOutput("home_sample_search_ui"),
          shiny::actionButton("home_open", "Open in Browser", class="btn-primary"))))))
}

rc_home_server <- function(input, output, session, all_exp, metadata, rv, df, libs) {
  catalog <- rc_home_catalog(all_exp, metadata)
  order <- c("Mammals","Other vertebrates","Invertebrates","Plants","Fungi","Protists","Bacteria","Archaea","Viruses","Unclassified")
  catalog <- catalog[order(match(catalog$clade,order),-catalog$libraries,catalog$name),,drop=FALSE]
  session$onFlushed(function() session$sendCustomMessage("rc-home-wheel", catalog), once=TRUE)
  shiny::observeEvent(input$home_wheel_species, {
    x <- as.character(input$home_wheel_species)
    if (length(x)==1L && x %in% catalog$name) shiny::updateSelectizeInput(session,"home_species",selected=x)
  }, ignoreInit=TRUE)
  shiny::observeEvent(input$home_species, {
    e <- all_exp[rc_home_host(organism)==input$home_species]
    groups <- rc_home_groups(e$name)
    ord <- order(match(groups,c("Modalities","Merged tracks","Collections","Single experiments")),e$name)
    choices <- stats::setNames(e$name[ord],paste(groups[ord],e$name[ord],sep=" · "))
    selected <- if (length(choices)) unname(choices[1]) else character()
    shiny::updateSelectizeInput(session,"home_experiment",choices=choices,selected=selected,server=TRUE)
  }, ignoreNULL=FALSE)
  output$home_species_info <- shiny::renderUI({
    x <- catalog[match(input$home_species,catalog$name),,drop=FALSE]; shiny::req(nrow(x),!is.na(x$name))
    shiny::div(class="rc-home-stats",sprintf("%s · %s Ribo-seq %s · %s %s",x$clade,
      format(x$libraries,big.mark=","),if (x$libraries == 1) "library" else "libraries",
      format(x$studies,big.mark=","),if (!is.na(x$studies) && x$studies == 1) "study" else "studies"))
  })
  output$home_sample_search_ui <- shiny::renderUI({
    shiny::req(input$home_experiment)
    if (!grepl("^all_samples-",input$home_experiment)) return(NULL)
    shiny::tagList(shiny::tags$hr(),shiny::h4("Find samples in this collection"),
      shiny::textInput("home_sample_query","Search sample metadata",placeholder="e.g. EIF4"),
      shiny::helpText("Searches collection annotations and supplied metadata. Select the samples to open; leave empty to use the default."),
      shiny::selectizeInput("home_samples","Matching samples",choices=character(),multiple=TRUE,width="100%",options=list(maxOptions=100)),
      shiny::textOutput("home_sample_count"))
  })
  sample_data <- shiny::eventReactive(input$home_experiment, {
    shiny::req(input$home_experiment,grepl("^all_samples-",input$home_experiment))
    tryCatch({ e <- rc_read_experiment(input$home_experiment,validate=FALSE,in.dir=ORFik::config()["exp"])
      base <- data.table::as.data.table(as.data.frame(e)); base[, Run := ORFik::runIDs(e)]
      base <- base[, intersect(c("Run","libtype","stage","condition","fraction","rep"),names(base)),with=FALSE]
      extra <- if (!is.null(metadata) && "Run" %in% names(metadata)) metadata[Run %in% base$Run] else NULL
      data.table::rbindlist(list(base,extra),fill=TRUE) },error=function(e) {warning(conditionMessage(e));data.table::data.table()})
  },ignoreNULL=TRUE)
  matches <- shiny::reactive({
    q <- trimws(input$home_sample_query %||% ""); shiny::req(nchar(q)>=2)
    m <- sample_data(); shiny::req(nrow(m))
    rc_home_search(m,q)
  })
  shiny::observeEvent(list(input$home_sample_query,input$home_experiment), {
    hits <- tryCatch(matches(),shiny.silent.error=function(e) character()); shiny::updateSelectizeInput(session,"home_samples",choices=stats::setNames(hits,hits),selected=character(),server=TRUE)
  },ignoreInit=TRUE)
  output$home_sample_count <- shiny::renderText({ q<-trimws(input$home_sample_query %||% ""); if(nchar(q)<2) return("Enter at least two characters."); sprintf("%s matching samples",format(length(matches()),big.mark=",")) })
  pending <- shiny::reactiveVal(NULL)
  shiny::observeEvent(input$home_open, {
    exp <- input$home_experiment; shiny::req(length(exp)==1L,exp %in% all_exp$name)
    ids <- if (grepl("^all_samples-",exp)) intersect(input$home_samples %||% character(),tryCatch(matches(),shiny.silent.error=function(e) character())) else character()
    pending(list(exp=exp,ids=ids)); rv$genome <- "ALL"; rv$exp <- exp
    shiny::updateNavbarPage(session,"navbarID",selected="browser")
  },ignoreInit=TRUE)
  shiny::observe({
    p <- pending(); shiny::req(!is.null(p)); current <- df(); shiny::req(identical(as.character(ORFik::name(current)),p$exp))
    available <- libs(); ids <- ORFik::runIDs(current); shiny::req(length(ids)==length(available))
    if (length(p$ids)) {
      chosen <- available[ids %in% p$ids]
      if (length(chosen)) session$onFlushed(function() shiny::updateSelectizeInput(session,"browser-library",choices=available,selected=chosen,server=TRUE),once=TRUE)
    }
    pending(NULL)
  },priority=-100)
}
