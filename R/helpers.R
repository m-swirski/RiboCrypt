fetchJS <- function(filename) {
    with_ext <- paste(filename, "js", sep = ".")
    path_to_file <- here("js", with_ext)
    lines <- readLines(path_to_file)
    paste(lines, sep = "", collapse = "")
}