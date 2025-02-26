lineDeSimplify <- function(plot) {
  plot <- plotly_build(plot)
  for (i in 1:length(plot$x$data)) plot$x$data[[i]]$line$simplify <- FALSE
  return(plot)
}

automateTicks <- function(plot) {
  plot %>% ggplotly(dynamicTicks = TRUE, tooltip=c("position")) %>%
    plotly::layout(yaxis=list(autorange = FALSE),xaxis=list(autorange=FALSE))
}

automateTicksLetters <- function(plot) {
  suppressWarnings(plot %>% ggplotlyHover(dynamicTicks = TRUE) %>%
                     plotly::layout(yaxis=list(autorange = FALSE),xaxis=list(autorange=FALSE)) %>%
                     toWebGL())
}


automateTicksGMP <- function(plot) {
  plot %>% ggplotly(dynamicTicks = TRUE, tooltip = "gene_names") %>%
    plotly::layout(yaxis=list(autorange = FALSE), xaxis=list(autorange=FALSE)) %>%
    style(hoverinfo = "text")
}
#'
#' @rawNamespace import(plotly, except = c(config, last_plot))
#' @keywords internal
automateTicksRNA <- function(plot, as_plotly = TRUE) {
  if (!as_plotly) return(plot)
  plot %>% ggplotly(dynamicTicks = TRUE) %>%
    plotly::layout(yaxis=list(autorange = FALSE, nticks=3), xaxis=list(autorange=FALSE))
}

automateTicksX <- function(plot) {
  plot %>% ggplotly(dynamicTicks = TRUE) %>%
    plotly::layout(xaxis=list(autorange = FALSE),
                   yaxis=list(autorange=FALSE,
                              ticktext=list("2","1","0"),
                              tickvals=list(-0.5,0.5,1.5),
                              tickmode = "array")
  ) %>% style(hoverinfo = "none")

}

#' Call ggplotly with hoveron defined
#' @param x a a ggplot argument
#' @param ... additional arguments for ggplotly
#' @return a ggplotly object
#' @keywords internal
ggplotlyHover <- function(x, ...) {
  gg <- plotly::ggplotly(x, ...)
  gg$x$data <- lapply(gg$x$data, function(x) {
    x$hoveron <- NULL
    x
  })
  return(gg)
}

addToImageButtonOptions <- function(multiomics_plot, filename, width, height,
                                    format = "svg") {
  multiomics_plot %>% plotly::config(
    toImageButtonOptions = list(
      format = format,
      filename = filename,
      width = width,
      height = height))
}

