fetchJS <- function(filename) {
  with_ext <- paste(filename, "js", sep = ".")
  path_to_file <- here::here("js", with_ext)
  lines <- readLines(path_to_file)
  paste(lines, sep = "", collapse = "")
}

get_jsdata <- function(ref_seq, display_range){
  seq <- extractTranscriptSeqs(findFa(ref_seq), display_range) %>% as.character %>% strsplit("") %>% .[[1]]
  js_data <- list(frame0 = seq[((seq_along(seq) - 1) %% 3 == 0)],
                  frame1 = seq[((seq_along(seq) - 1) %% 3 == 1)],
                  frame2 = seq[((seq_along(seq) - 1) %% 3 == 2)],
                  colors = c("#F8766D","#00BA38","#619CFF"))
  js_data
}

