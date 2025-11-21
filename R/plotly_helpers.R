lineDeSimplify <- function(plot) {
  plot <- plotly_build(plot)
  for (i in 1:length(plot$x$data)) {
    if (!is.null(plot$x$data[[i]]$line)) plot$x$data[[i]]$line$simplify <- FALSE
  }
  return(plot)
}

automateTicks <- function(plot) {
  plot %>% ggplotly(dynamicTicks = TRUE, tooltip=c("position")) %>%
    plotly::layout(yaxis=list(autorange = FALSE), xaxis=list(autorange=FALSE))
}

automateTicksLetters <- function(plot) {
  suppressWarnings(
    plot %>% ggplotlyHover(dynamicTicks = TRUE) %>%
                     plotly::layout(yaxis=list(autorange = FALSE, fixedrange = TRUE),
                                    xaxis=list(autorange=FALSE)) %>%
                     toWebGL())
}

automateTicksGMP <- function(plot) {
  suppressWarnings(
    plot %>% ggplotly(dynamicTicks = TRUE, tooltip = "gene_names") %>%
      plotly::layout(yaxis=list(autorange = FALSE), xaxis=list(autorange=FALSE)) %>%
      style(hoverinfo = "text"))
}

automateTicksAA <- function(plot) {
  suppressWarnings(
    plot %>% ggplotlyHover(dynamicTicks = TRUE) %>%
      plotly::layout(yaxis=list(autorange = FALSE, fixedrange = TRUE),
                     xaxis=list(autorange=FALSE)) %>%
      toWebGL())
}

nice_ticks <- function(y_max, n_ticks = 3) {
  if (y_max <= 0) return(numeric(0))

  # Ask pretty() for one extra tick so we have room to drop 0
  raw <- pretty(c(0, y_max), n = n_ticks + 1)

  # Drop 0 and negatives
  raw_pos <- raw[raw > 0]

  # Keep at most n_ticks, biased toward the top
  if (length(raw_pos) > n_ticks) {
    raw_pos_even <- raw_pos[seq_along(raw_pos) %% 2 == 0]
    raw_pos <- c(head(raw_pos_even, length(raw_pos_even) - 1), tail(raw_pos, 1))
  }

  raw_pos
}


#'
#' @rawNamespace import(plotly, except = c(config, last_plot))
#' @keywords internal
automateTicksRNA <- function(plot, as_plotly = TRUE, scale_ticks = TRUE) {
  if (!as_plotly) return(plot)
  if (scale_ticks) {
    y_max <- max(plot$data$count)
    y_nticks <- ifelse(y_max > 3,
                       ifelse(y_max > 10, 3, 2), 1)
    tick_vals <- nice_ticks(y_max, n_ticks = y_nticks)
    plot %>% ggplotly() %>%
      plotly::layout(yaxis=list(range = c(0, max(tick_vals)),
                                tickvals = tick_vals,
                                ticktext = tick_vals),
                     xaxis=list(autorange=FALSE))
  } else {
    plot %>% ggplotly() %>%  plotly::layout(xaxis=list(autorange=FALSE))
  }
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
  gg <- plotly::ggplotly(x, ..., tooltip = "text")
  gg$x$data <- lapply(gg$x$data, function(x) {
    x$hoveron <- NULL
    x
  })
  return(gg)
}

addToImageButtonOptions <- function(multiomics_plot, filename, width, height,
                                    format = "svg", modeBarButtonsToRemove = c("lasso2d", "select2d", "zoomIn2d", "zoomOut2d",
                                                                               "hoverCompareCartesian", "hoverClosestCartesian")) {
  # TODO: Add a proper button
  ribocrypt_svg_path <- "M24.92 12.183c0-1.586-.604-2.864-1.585-3.83.172-.547.398-1.763-.229-3.321 0 0-1.114-.348-3.628 1.315a12.695 12.695 0 0 0-3.081-.366c-1.154 0-2.322.143-3.409.44-2.596-1.747-3.74-1.391-3.74-1.391-.748 1.847-.287 3.215-.145 3.554-.883.936-1.414 2.133-1.414 3.594 0 1.111.128 2.099.44 2.964l.325.732c.879 1.614 2.606 2.655 5.677 2.983-.434.289-.885.779-1.062 1.612-.594.28-2.475.966-3.603-.944 0 0-.633-1.148-1.842-1.235 0 0-1.174-.017-.08.722 0 0 .782.367 1.326 1.738 0 0 .705 2.342 4.114 1.593v2.417s-.076.857-.867 1.143c0 0-.469.312.034.497 0 0 2.205.174 2.205-1.604v-2.643s-.09-1.047.429-1.404v4.332s-.032 1.031-.576 1.421c0 0-.362.646.433.468 0 0 1.517-.211 1.584-1.967l.035-4.383h.363l.033 4.383c.076 1.748 1.59 1.967 1.59 1.967.793.179.429-.468.429-.468-.54-.389-.579-1.421-.579-1.421v-4.297c.52.402.436 1.369.436 1.369v2.643c0 1.777 2.2 1.604 2.2 1.604.505-.186.036-.498.036-.498-.793-.286-.867-1.143-.867-1.143v-3.461c0-1.346-.574-2.056-1.137-2.435 3.277-.318 4.845-1.368 5.572-2.99-.015.027.26-.726.26-.726.25-.859.325-1.855.325-2.963h-.002z"

  ribocrypt <- list(
    name = "RiboCrypt",
    icon = list(
      path = ribocrypt_svg_path,  # Using absolute URL format for Plotly
      transform = 'matrix(1 0 0 1 -2 -2) scale(0.7)'  # Optional scaling and positioning
    ),
    click = htmlwidgets::JS(
      "function(gd) {
       var txt = {x: [1], y: [1], text: 'RiboCrypt', mode: 'text'};
       Plotly.addTraces(gd, txt);
    }"
    )
  )

  multiomics_plot %>% plotly::config(
    toImageButtonOptions = list(
      format = format,
      filename = filename,
      width = width,
      height = height
    ),
    displaylogo = FALSE,
    modeBarButtonsToAdd = list(ribocrypt),
    modeBarButtonsToRemove = modeBarButtonsToRemove
  )
}

