fetchJS <- function(script_name) {
  script = system.file("js", script_name, package =
                                                "RiboCrypt")
  lines <- readLines(script)
  paste(lines, sep = "", collapse = "")
}

fetch_JS_seq <- function(target_seq, nplots, distance = 50, display_dist) {
  nt_yaxis <- paste("y", nplots + 1, sep = "")
  x1 <- seq(1,display_dist,3)
  x2 <- seq(2,display_dist,3)
  x3 <- seq(3,display_dist,3)
  rendered_seq <- strsplit(as.character(target_seq),"")[[1]]
  js_data <- list(list(x = x1,
                       text =  rendered_seq[x1],
                       y = rep(0.5, length(x1)),
                       xaxis = "x",
                       yaxis = nt_yaxis,
                       distance = distance,
                       color = "#F8766D"
                       ),
                  list(x = x2,
                       text = rendered_seq[x2],
                       y = rep(0.5,length(x2)),
                       xaxis = "x",
                       yaxis = nt_yaxis,
                       distance = distance,
                       color = "#00BA38"
                       ),
                  list(x = x3,
                       text = rendered_seq[x2],
                       y = rep(0.5, length(x3)),
                       xaxis = "x",
                       yaxis = nt_yaxis,
                       distance = distance,
                       color = "#619CFF"
                       ))

  js_data
}
