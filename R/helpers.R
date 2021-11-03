fetchJS <- function(path_to_file) {
    lines <- readLines(path_to_file)
    paste(lines, sep = "", collapse = "")
}
