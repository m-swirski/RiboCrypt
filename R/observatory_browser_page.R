observatory_browser_ui <- function(id) {
  ns <- shiny::NS(id)

  shiny::tabPanel(
    "Browse",
    shiny::fluidRow(
      shiny::column(
        2,
        shiny::selectizeInput(ns("gene_input"),
          "Gene",
          choices = list()
        )
      ),
      shiny::column(
        2,
        shiny::selectizeInput(
          ns("tx_input"),
          "Transcript",
          choices = list()
        )
      ),
      shiny::column(1, plot_button(ns("go")))
    )
  )
}

observatory_browser_server <- function(
  id,
  meta_experiment_df,
  library_selections
) {
  shiny::moduleServer(id, function(input, output, session) {
    # Derive gene <-> transcript mapping table from the current experiment.
    # Returns a data.table with columns: value (tx id), label (gene symbol).
    gene_name_list <- shiny::reactive({
      get_gene_name_categories(meta_experiment_df())
    }) |> shiny::bindEvent(meta_experiment_df(), ignoreNULL = TRUE)

    # When the experiment changes, repopulate the gene dropdown with the new
    # gene symbols and pre-select the first one.
    shiny::observe({
      gene_names <- unique(gene_name_list()[, 2][[1]])
      shiny::updateSelectizeInput(
        session,
        "gene_input",
        choices = gene_names,
        selected = gene_names[1],
        server = TRUE
      )
    }) |> shiny::bindEvent(gene_name_list(), ignoreNULL = TRUE)

    # When the selected gene changes, repopulate the transcript dropdown with
    # only the isoforms that belong to the chosen gene.
    shiny::observe({
      shiny::req(input$gene_input != "")
      isoforms <- gene_name_list()[label == input$gene_input, value]
      shiny::updateSelectizeInput(
        session,
        "tx_input",
        choices = isoforms,
        selected = isoforms[1],
        server = TRUE
      )
    }) |> shiny::bindEvent(
      input$gene_input,
      ignoreNULL = TRUE, ignoreInit = TRUE
    )
  })
}
