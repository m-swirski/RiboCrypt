lineDeSimplify <- function(plot) {
  plot <- plotly_build(plot)
  for (i in 1:length(plot$x$data)) plot$x$data[[i]]$line$simplify <- FALSE
  return(plot)
}

automateTicks <- function(plot) {
  plot %>% ggplotly(dynamicTicks = TRUE) %>% plotly::layout(yaxis=list(autorange = FALSE),xaxis=list(autorange=FALSE))
}
#'
#' @rawNamespace import(plotly, except = c(config, last_plot))
#' @keywords internal
automateTicksRNA <- function(plot) {
  plot %>% ggplotly(dynamicTicks = TRUE) %>% plotly::layout(yaxis=list(autorange = FALSE,nticks=3),xaxis=list(autorange=FALSE))
}

automateTicksX <- function(plot) {
  plot %>% ggplotly(dynamicTicks = TRUE) %>% plotly::layout(xaxis=list(autorange = FALSE),
                                                            yaxis=list(autorange=FALSE,
                                                                       ticktext=list("2","1","0"),
                                                                       tickvals=list(-0.5,0.5,1.5),
                                                                       tickmode = "array")
  )
}
