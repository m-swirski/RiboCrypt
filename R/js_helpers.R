fetchJS <- function(script_name) {
  script = system.file("js", script_name, package =
                                                "RiboCrypt")
  lines <- readLines(script)
  paste(lines, sep = "", collapse = "")
}

fetch_JS_seq <- function(target_seq, nplots, distance = 50, display_dist) {
  fr_colors <- c("#F8766D","#00BA38", "#619CFF")
  nt_yaxis <- paste("y", nplots + 1, sep = "")
  aa_yaxis <- paste("y", nplots + 3, sep = "")
  nts <- sapply(1:3, function(x) seq(x, display_dist, 3))
  rendered_seq <- strsplit(as.character(target_seq),"")[[1]]
  aas <- lapply(1:3, function(x) suppressWarnings(strsplit(as.character(translate(target_seq[[1]][x:display_dist])), "")[[1]]))
  nts_js_data <- lapply(1:3, function(fr) list(x = nts[[fr]],
                                               text = rendered_seq[nts[[fr]]],
                                               y = rep(0.5, length(nts[[fr]])),
                                               xaxis = "x",
                                               yaxis = nt_yaxis,
                                               distance = distance,
                                               color = fr_colors[fr]))
  aa_js_data <- lapply(1:3, function(fr) list(x = nts[[fr]][seq_along(aas[[fr]])],
                                              text = aas[[fr]],
                                              y = rep(2.4 - fr, length(aas[[fr]])),
                                              xaxis = "x",
                                              yaxis = aa_yaxis,
                                              distance = distance * 2,
                                              color = "grey45"))

  return(c(nts_js_data, aa_js_data))
}
