# Narrow wrappers around ORFik internals used by RiboCrypt.
orfik_name_decider <- function(...) {
  getFromNamespace("name_decider", "ORFik")(...)
}

orfik_remove_file_ext <- function(...) {
  getFromNamespace("remove.file_ext", "ORFik")(...)
}
