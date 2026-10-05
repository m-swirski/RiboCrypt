#' Library membership for currently visible groups, independent of collapse.
#' @noRd
megabrowser_group_membership <- function(table, meta, selected = character()) {
  groups <- megabrowser_display_groups(meta, colnames(table))
  if (length(selected)) groups <- groups[intersect(names(groups), selected)]
  rows <- lapply(names(groups), function(group) {
    members <- copy(meta[match(colnames(table)[groups[[group]]], Run)])
    members[, Group := group]
    members[, setdiff(names(members), c("index", "order", "cluster")), with = FALSE]
  })
  rbindlist(rows, use.names = TRUE)
}

#' Summarize metadata without treating a mode as an enrichment result.
#' @noRd
megabrowser_metadata_summary <- function(values) {
  missing <- is.na(values)
  numeric <- is.numeric(values)
  if (!numeric) missing <- missing | as.character(values) == ""
  present <- values[!missing]
  value <- if (length(present)) megabrowser_group_value(present) else NA
  list(Value = as.character(value), `Agreement (%)` = if (numeric || !length(present)) NA_real_ else
         100 * sum(as.character(present) == value) / length(values), Missing = sum(missing),
       Minimum = if (numeric && length(present)) min(present) else NA_real_,
       Maximum = if (numeric && length(present)) max(present) else NA_real_)
}

#' One row per group and metadata field, always using original library values.
#' @noRd
megabrowser_group_summary <- function(membership) {
  fields <- setdiff(names(membership), c("Group", "Run"))
  rbindlist(lapply(unique(membership$Group), function(group) {
    members <- membership[Group == group]
    rbindlist(lapply(fields, function(field) {
      cbind(data.table(Group = group, Libraries = nrow(members), Field = field),
            as.data.table(megabrowser_metadata_summary(members[[field]])))
    }))
  }))
}

#' Downloads and a lazy group-inspection table.
#' @noRd
megabrowser_group_outputs <- function(output, table, grouped, selected) {
  membership <- reactive(megabrowser_group_membership(table()$table, grouped()$meta, selected()))
  summary <- reactive(megabrowser_group_summary(membership()))
  output$group_summary <- DT::renderDT(DT::datatable(summary(), rownames = FALSE,
    options = list(pageLength = 15, scrollX = TRUE)) %>%
      DT::formatRound(c("Agreement (%)", "Minimum", "Maximum"), 2), server = TRUE)
  output$download_groups <- megabrowser_csv_download(summary, "group-summary")
  output$download_membership <- megabrowser_csv_download(membership, "library-membership")
}

#' CSV values retain their precision; rounding is for the table only.
#' @noRd
megabrowser_csv_download <- function(data, name) {
  downloadHandler(filename = function() paste0("megabrowser-", name, ".csv"),
    content = function(file) data.table::fwrite(data(), file, na = ""), contentType = "text/csv")
}
