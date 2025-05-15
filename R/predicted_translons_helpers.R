# Function to load data
load_data <- function(species) {
  df <- read.experiment(species, validate = FALSE)
  table_path <- file.path(refFolder(df),
                          "predicted_translons",
                          "predicted_translons_with_sequence.fst")
  if (file.exists(table_path)) {
    translon_table <- fst::read_fst(table_path, as.data.table = TRUE)
    setattr(translon_table, "exp", species)
  } else {
    NULL
  }

  reactiveValues(translon_table = translon_table, df = df)
}

# Function to load data
load_data_umap <- function(species) {
  df <- read.experiment(species, validate = FALSE)
  dir <- file.path(refFolder(df), "UMAP")
  table_path <- file.path(dir, "UMAP_by_gene_counts.fst")
  if (file.exists(table_path)) {
    dt_umap <- fst::read_fst(table_path, as.data.table = TRUE)
    setattr(dt_umap, "exp", species)
  } else {
    NULL
  }

  reactiveValues(dt_umap = dt_umap, df = df)
}

# Generalized function to handle download trigger buttons.
# 'format' should be either "csv" or "xlsx".
# 'trigger_input' is the name of the visible actionButton (e.g. "trigger_download_csv").
# 'download_button' is the id (as a string) of the hidden downloadButton (e.g. "download_csv").
handle_download_trigger <- function(input, output, current_format, trigger_input, download_button, md, session) {
  with(rlang::caller_env(), {
    observeEvent(input[[trigger_input]], {
      if (is.null(md()) || isolate(input$dff) != name(md()$df)) md(load_data(isolate(input$dff)))
      req(md()$df)
      # For Excel, check that a table is available (otherwise abort)
      if (is.null(md()$translon_table)) {
        showNotification("No translon predictions for this dataset!", type = "warning")
        req(md()$translon_table)
      }
      # Compute the filename for this dataset
      filename <- generate_filename(md()$df, current_format, FALSE)
      current_downloads <- downloaded_files()
      if (filename %in% current_downloads) {
        showNotification(
          paste("You have already downloaded the", ifelse(current_format == "csv", "CSV", "Excel"),
                "for this dataset!"),
          type = "warning"
        )
      } else {
        downloaded_files(c(current_downloads, filename))
        message("Clicking hidden button for ", current_format)
        shinyjs::click(download_button)
      }
    })
  })
}

# Generalized function for download handlers.
# 'write_fun' is a function that writes the table to a file.
make_download_handler <- function(format, write_fun, md) {
  downloadHandler(
    filename = function() {
      generate_filename(md()$df, format)
    },
    content = function(file) {
      showNotification(paste("Preparing", format, "file.."), type = "message")
      write_fun(file)
    }
  )
}

generate_filename <- function(df, format, show_message = TRUE) {
  file <- gsub(" ", "_", paste0(organism(df), "_predicted_translons_RiboCrypt.", format))
  if (show_message) message("Generated filename: ", file)
  file
}

# Generalized function to render the DT table.
render_translon_datatable <- function(data, session, add_links = TRUE) {
  # Add URL
  if (add_links) {
    host <- getHostFromURL(session) #"https://ribocrypt.org"
    exp <- attr(data, "exp")
    # exp <- "all_merged-Homo_sapiens_modalities"
    symbols <- if(!is.null(data$external_gene_name)) {data$external_gene_name} else NULL
    gene_ids <- if(!is.null(data$ensembl_gene_id)) {data$ensembl_gene_id} else data$GENE
    tx_ids <- if(!is.null(data$ensembl_tx_name)) {data$ensembl_tx_name} else data$TX
    urls <- make_rc_url(symbol = symbols, gene_id = gene_ids, tx_id = tx_ids,
                        exp = exp,
                        libraries = NULL, leader_extension = 2000, trailer_extension = 0,
                        viewMode = FALSE, other_tx = FALSE,
                        plot_on_start = TRUE, frames_type = "columns", kmer=1,
                        add_translons = TRUE, zoom_range = data$coordinates,
                        host = host)
    data <- cbind(link = paste0('<a href=', urls, ' target=\"_blank\">',
                                "Translon_", seq(length(urls)), '</a>'),
                  data)
  }

  datatable(data,
            extensions = 'Buttons',
            escape = FALSE,
            filter = "top",
            options = list(
              dom = 'Bfrtip',
              buttons = list(
                list(
                  extend = "csv",
                  text = "Download current page (CSV)",
                  filename = "current",
                  exportOptions = list(
                    modifier = list(page = "current")
                  )
                ),
                list(
                  extend = "excel",
                  text = "Download current page (Excel)",
                  filename = "current",
                  exportOptions = list(
                    modifier = list(page = "current")
                  )
                )
              )
            )
  )
}
