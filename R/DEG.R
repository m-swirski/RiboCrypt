DEG_plot <- function(dt) {
  dt <- highlight_key(dt, ~id, "Select a transcript")

  gg <- ggplot(dt, aes(x = rna, y = rfp, color = Regulation, frame = variable)) +
    geom_vline(aes(xintercept = 0), color = "red")+
    geom_hline(aes(yintercept = 0), color = "red") +
    geom_abline(slope = 1, intercept = 0, color = "grey", linetype = 2)+
    geom_point() +
    xlim(c(-1,1))+
    ylim(c(-1,1))

  select <- highlight(
    ggplotly(gg),
    on = "plotly_selected",
    selectize = TRUE, persistent = TRUE
  )
  return(select)
}
