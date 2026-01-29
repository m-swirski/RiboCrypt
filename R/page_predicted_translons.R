predicted_translons_ui <- function(id, all_exp_translons, label = "predicted_translons") {
  ns <- NS(id)
  tabPanel(
    title = "Predicted Translons", icon = icon("rectangle-list"),
    h2("Predicted Translons Overview"),
    # Include shinyjs so we can trigger hidden buttons
    shinyjs::useShinyjs(),
    sidebarLayout(
      sidebarPanel(
        tags$style(HTML("
        table.dataTable td.dt-id {
          color: #007BFF;       /* Bootstrap link blue */
          cursor: pointer;      /* hand cursor on hover */
          text-decoration: underline;
        }
        table.dataTable td.dt-id:hover {
          color: #0056b3;       /* darker on hover */
        }
      ")),
        experiment_input_select(all_exp_translons$name, ns),
        actionButton(ns("go"), "Search", icon = icon("magnifying-glass")),
        hr(),
        tags$b("Download Full Table:"),
        div(
          # Visible download trigger buttons
          actionButton(ns("trigger_download_csv"), "Download CSV",
                       icon = icon("file-csv"), class = "btn btn-primary"),
          actionButton(ns("trigger_download_excel"), "Download Excel",
                       icon = icon("file-excel"), class = "btn btn-success"),
          # Hidden download buttons that use downloadHandler
          downloadButton(ns("download_csv"), label = NULL, style = "visibility: hidden;"),
          downloadButton(ns("download_excel"), label = NULL, style = "visibility: hidden;"),
          div(style = "display:none;",
            checkboxInput(ns("useCustomRegions"), "Protein structures", TRUE)
          ),
          div(style = "display:none;",
            textInput(ns("selectedRegion"), NULL, value = "")
          ),
          style = "display: flex; gap: 10px; margin-top: 10px;"
        )
      ),
      mainPanel(
        DT::DTOutput(ns("translon_table")) %>% shinycssloaders::withSpinner(color = "#0dc5c1"),
        uiOutput(ns("proteinStruct"))
      )
    )
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
        df <- read.experiment(attr(md()$translon_table, "exp"), validate = FALSE)
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
