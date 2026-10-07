#' Experiments supported by the Predicted Translons page.
#' @noRd
predicted_translons_experiments <- function(all_exp) {
  all_exp[grep("all_merged", name), ][libtypes == "RFP", ][
    grep("Escherichia_coli", name, invert = TRUE), ]
}

predicted_translons_ui <- function(id, all_exp_translons, label = "predicted_translons") {
  ns <- NS(id)
  tabPanel(
    title = "Translons", value = "Predicted Translons", icon = icon("rectangle-list"),
    h2("Predicted Translons Overview"),
    # Include shinyjs so we can trigger hidden buttons
    shinyjs::useShinyjs(),
    tags$style(HTML("
        .rc-translon-toolbar { display: flex; flex-wrap: wrap; align-items: center; gap: 16px 48px; padding: 8px 0 12px; margin-bottom: 12px; }
        .rc-translon-study-group { flex: 1 1 360px; min-width: 0; max-width: 560px; }
        .rc-translon-study { display: flex; align-items: center; gap: 8px; min-width: 0; }
        .rc-translon-study .shiny-input-container { flex: 1; min-width: 0; width: auto !important; margin: 0 !important; }
        .rc-translon-toolbar .selectize-control { margin-bottom: 0; }
        .rc-translon-simplified .shiny-input-container { width: auto; margin: 0 !important; }
        .rc-translon-toolbar .checkbox { margin: 0; }
        .rc-translon-downloads { display: grid; gap: 8px; }
        .rc-translon-downloads .rc-translon-excel { background-color: #217346 !important; border-color: #217346 !important; color: white !important; }
        .rc-translon-downloads .rc-translon-excel:hover,
        .rc-translon-downloads .rc-translon-excel:focus { background-color: #185c37 !important; border-color: #185c37 !important; }
        table.dataTable td.dt-id {
          color: #007BFF;       /* Bootstrap link blue */
          cursor: pointer;      /* hand cursor on hover */
          text-decoration: underline;
        }
        table.dataTable td.dt-id:hover {
          color: #0056b3;       /* darker on hover */
        }
      ")),
    div(class = "rc-translon-toolbar",
      div(class = "rc-translon-study-group",
        tags$label("Study", `for` = ns("dff")),
        div(class = "rc-translon-study",
          experiment_input_select(all_exp_translons$name, ns, label = NULL),
          actionButton(ns("go"), "Search", icon = icon("magnifying-glass")))),
      div(class = "rc-translon-simplified", checkboxInput(ns("simplified"), "simplified", FALSE)),
      div(class = "rc-translon-downloads",
        actionButton(ns("trigger_download_csv"), "Download full CSV",
                     icon = icon("file-csv"), class = "btn btn-primary"),
        actionButton(ns("trigger_download_excel"), "Download full Excel",
                     icon = icon("file-excel"), class = "btn btn-success rc-translon-excel"))
    ),
    div(style = "display:none;",
      downloadButton(ns("download_csv"), label = NULL),
      downloadButton(ns("download_excel"), label = NULL),
      checkboxInput(ns("useCustomRegions"), "Protein structures", TRUE),
      textInput(ns("selectedRegion"), NULL, value = "")
    ),
    DT::DTOutput(ns("translon_table")) %>% shinycssloaders::withSpinner(color = "#0dc5c1"),
    uiOutput(ns("proteinStruct"))
  )
}


predicted_translons_server <- function(id, all_exp, browser_options) {
  moduleServer(
    id,
    function(input, output, session) {
      # Track if "Plot" has been clicked
      plot_triggered <- reactiveVal(FALSE)
      download_trigger <- reactiveVal(NULL)
      # Reactive vector to store downloaded filenames (for full dataset downloads)
      downloaded_files <- reactiveVal(character(0))
      # Reactive to store data (loads when needed)
      md <- reactiveVal(NULL)

      # Set human to auto, clean this code later

      # Trigger data loading when "Plot" is clicked.
      # Also reset the downloaded files vector.
      observeEvent(input$go, {
        req(isTruthy(input$dff))
        md(load_data(isolate(input$dff)))
        plot_triggered(TRUE)
      })


      # Render DT Table ONLY if "Plot" was clicked
      output$translon_table <- DT::renderDT({
        req(plot_triggered())
        render_translon_datatable(md()$translon_table, session)
      }, server = TRUE)
      observeEvent(download_trigger(), {
        req(download_trigger())  # Ensure a value is set
        type <- download_trigger()
        trigger_input <- paste0("trigger_download_", type)
        message("Firing button: ", trigger_input)

        # Fire the actual button event
        Sys.sleep(0.6)
        shinyjs::click(trigger_input)

        # Reset download trigger after firing
        download_trigger(NULL)
      }, priority = -300)

      observeEvent(input$translon_id_click, {
        id  <- input$translon_id_click$id
        req(id != "")
        row <- input$translon_id_click$row
        df <- rc_read_experiment(attr(md()$translon_table, "exp"), validate = FALSE)
        pep_dir <- file.path(refFolder(df), "protein_structure_predictions")
        path <- pep_id_to_path(id, pep_dir)
        path <- ifelse(path == "", "No file found", path)
        showNotification(paste0("Clicked Protein id: ", id," (", path, ")"))
        if (path != "No file found") {
          updateTextInput(inputId = "selectedRegion", value = path)
          shinyjs::runjs("window.scrollTo(0, document.body.scrollHeight);")
        }
      })
      gene_name_list <- reactiveVal(data.table())
      module_protein(input, output, gene_name_list, session)

      translon_specific_url_checker()

      selected_exp <- browser_options["default_experiment_translon"]
      experiment_update_select(NULL, all_exp, all_exp$name, selected_exp)

      # Use the helper for both CSV and Excel
      for (format in c("csv", "xlsx")) {
        local({
          current_format <- format
          type <- ifelse(current_format == "xlsx", "excel", current_format)
          trigger_input <- paste0("trigger_download_", type)
          download_button <- paste0("download_", type)
          handle_download_trigger(input, output, current_format, trigger_input, download_button, md, session)
          output[[download_button]] <- make_download_handler(current_format, function(file) {
            if (current_format == "xlsx") {
              write_xlsx(md()$translon_table, file)
            } else if (current_format == "csv") {
              fwrite(md()$translon_table, file)
            } else {
              stop("Invalid format for translon download!")
            }
          }, md)
          outputOptions(output, download_button, suspendWhenHidden = FALSE)
        })
      }
    }
  )
}
