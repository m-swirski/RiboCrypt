library_selection_picker <- function(id) {
  ns <- shiny::NS(id)

  shiny::tagList(
    shiny::selectizeInput(
      ns("active_selection_id"),
      "Selection",
      choices = list()
    ),
    shiny::textInput(
      ns("active_selection_label"),
      "Selection name",
      value = ""
    )
  )
}

library_selection_reset_button <- function(id) {
  ns <- shiny::NS(id)
  shiny::actionButton(ns("reset_active_selection"), "Reset")
}

# reactive_plot_selection -
# a reactiveVal holding a vector of character vectors representing runIds
# reactive_dataplot_selection -
# a reactiveVal holding a vector of character vectors representing runIds
# reactive_dataplot_selection is always a subset of reactive_plot_selection
library_selection_server <- function(
  id,
  reactive_plot_selection_input,
  reactive_data_table_selection_input
) {
  shiny::moduleServer(id, function(input, output, session) {
    # label for a choice that will result in creating a new selection
    new_selection_choice <- "New selection..."

    # counter used to generate unique ids for selections
    counter <- shiny::reactiveVal(2)
    # reactive value to store multiple selections
    selection_store <- {
      plot_selections <- list()
      plot_selections[[as.character(1)]] <- NULL
      data_table_selections <- list()
      data_table_selections[[as.character(1)]] <- NULL
      labels <- list()
      labels[[as.character(1)]] <- ""
      list(
        index = shiny::reactiveVal(c(as.character(1))),
        plot_selections = shiny::reactiveVal(plot_selections),
        data_table_selections = shiny::reactiveVal(data_table_selections),
        labels = shiny::reactiveVal(labels)
      )
    }

    # reactive value to keep track of which selection is active at the moment
    active_selection_id <- shiny::reactiveVal(as.character(1))

    # reactive value to keep track of the current plot selection
    active_plot_selection <- shiny::reactive({
      shiny::req(!is.null(active_selection_id()) && active_selection_id() != "")
      selection_store$plot_selections()[[active_selection_id()]]
    }) |> shiny::bindEvent(
      selection_store$plot_selections(),
      active_selection_id()
    )

    # reactive value to keep track of the current data_table selection
    active_data_table_selection <- shiny::reactive({
      shiny::req(!is.null(active_selection_id()) && active_selection_id() != "")
      selection_store$data_table_selections()[[active_selection_id()]]
    }) |> shiny::bindEvent(
      selection_store$data_table_selections(),
      active_selection_id()
    )

    # Observer to keep input in sync with selection store
    shiny::observe({
      labels <- selection_store$labels()
      selection_ids <- selection_store$index()
      display_labels <- vapply(
        selection_ids,
        function(selection_id) {
          label <- labels[[selection_id]]
          if (is.null(label) || label == "" || label == selection_id) {
            selection_id
          } else {
            paste(selection_id, label, sep = " - ")
          }
        },
        character(1)
      )
      display_choices <- stats::setNames(selection_ids, display_labels)
      shiny::updateSelectizeInput(
        session,
        "active_selection_id",
        choices = c(display_choices, new_selection_choice),
        selected = active_selection_id()
      )
    }) |> shiny::bindEvent(selection_store$index(), selection_store$labels())

    # Keep selection name input in sync with active selection
    shiny::observe({
      label <- selection_store$labels()[[active_selection_id()]]
      shiny::updateTextInput(
        session,
        "active_selection_label",
        value = if (is.null(label)) "" else label
      )
    }) |> shiny::bindEvent(
      active_selection_id(),
      selection_store$labels()
    )

    # Store updated label for active selection
    shiny::observe({
      shiny::req(!is.null(active_selection_id()) && active_selection_id() != "")
      labels <- selection_store$labels()
      labels[[active_selection_id()]] <- input$active_selection_label
      selection_store$labels(labels)
    }) |> shiny::bindEvent(input$active_selection_label, ignoreInit = TRUE)

    set_active_label <- function(label, selection_id = active_selection_id()) {
      if (is.null(selection_id) || selection_id == "") return(invisible(NULL))
      if (is.null(label)) label <- ""
      labels <- selection_store$labels()
      labels[[selection_id]] <- label
      selection_store$labels(labels)
      if (selection_id == active_selection_id()) {
        shiny::updateTextInput(session, "active_selection_label", value = label)
      }
      invisible(NULL)
    }

    # Observer for creating a new selection
    shiny::observe({
      shiny::req(input$active_selection_id == new_selection_choice)
      new_selection_id <- counter()
      counter(new_selection_id + 1)

      # add new empty selection to the store
      index <- c(selection_store$index(), as.character(new_selection_id))

      plot_selections <- selection_store$plot_selections()
      plot_selections[[as.character(new_selection_id)]] <- NULL

      data_table_selections <- selection_store$data_table_selections()
      data_table_selections[[as.character(new_selection_id)]] <- NULL
      labels <- selection_store$labels()
      labels[[as.character(new_selection_id)]] <- ""

      # update the store
      selection_store$index(index)
      selection_store$plot_selections(plot_selections)
      selection_store$data_table_selections(data_table_selections)
      selection_store$labels(labels)

      # set new selection as active
      active_selection_id(as.character(new_selection_id))
    }) |> shiny::bindEvent(input$active_selection_id)

    # Observer for changing the active selection
    shiny::observe({
      shiny::req(
        !is.null(input$active_selection_id) &&
          input$active_selection_id != "" &&
          input$active_selection_id != new_selection_choice
      )
      active_selection_id(input$active_selection_id)
    }) |> shiny::bindEvent(input$active_selection_id)

    # Observer for reseting active selection
    shiny::observe({
      plot_selections <- selection_store$plot_selections()
      plot_selections[[active_selection_id()]] <-
        NULL
      selection_store$plot_selections(plot_selections)

      data_table_selections <- selection_store$data_table_selections()
      data_table_selections[[active_selection_id()]] <-
        NULL
      selection_store$data_table_selections(data_table_selections)
    }) |> shiny::bindEvent(
      input$reset_active_selection
    )

    # Observer to keep active plot_selection
    # in sync with incoming values of reactive_plot_selection
    shiny::observe({
      plot_selections <- selection_store$plot_selections()
      plot_selections[[active_selection_id()]] <-
        reactive_plot_selection_input()

      selection_store$plot_selections(plot_selections)
    }) |> shiny::bindEvent(reactive_plot_selection_input())

    # Observer to keep active_data_table_selection
    # in sync with incoming values of reactive_data_table_selection
    shiny::observe({
      data_table_selections <- selection_store$data_table_selections()
      data_table_selections[[active_selection_id()]] <-
        reactive_data_table_selection_input()
      selection_store$data_table_selections(data_table_selections)
    }) |> shiny::bindEvent(
      reactive_data_table_selection_input(),
      ignoreNULL = FALSE
    )

    list(
      active_selection_id = active_selection_id,
      active_plot_selection = active_plot_selection,
      active_data_table_selection = active_data_table_selection,
      set_active_label = set_active_label,
      all_selections = selection_store
    )
  })
}
