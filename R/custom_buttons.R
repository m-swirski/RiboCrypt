plot_button <- function(id, text = "Generate Plot", icon_type= icon("rocket"),
                        style = "width: 100%; font-size: 16px;") {
  actionButton(id, text, icon = icon_type,
               class = "btn-primary", style = style)
}
