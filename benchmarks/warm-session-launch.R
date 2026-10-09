# Source-loaded developer launcher; never evaluates commands after the first app.
devtools::load_all(".")
if (dir.exists(path.expand("~/Desktop/ORFik"))) devtools::load_all(path.expand("~/Desktop/ORFik"))
for (expr in parse(path.expand(Sys.getenv("RIBOCRYPT_RUN_SCRIPT", "~/Desktop/run_ribocrypt.R")))) {
  if (!is.call(expr)) next
  fun <- as.character(expr[[1]])
  if (identical(fun, "RiboCrypt_app")) {
    expr$init_tab_focus <- "browser"
    app <- eval(expr)
    break
  }
  if (identical(fun, "library") || any(fun == "load_all")) next
  eval(expr)
}
shiny::runApp(app, host = "127.0.0.1", port = as.integer(Sys.getenv("PORT", "7820")), launch.browser = FALSE)
