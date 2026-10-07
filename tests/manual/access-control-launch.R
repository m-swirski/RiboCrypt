# Test fixtures only: do not use these identities/secrets in production.
if (dir.exists(path.expand("~/Desktop/ORFik"))) devtools::load_all(path.expand("~/Desktop/ORFik"))
devtools::load_all(Sys.getenv("RIBOCRYPT_REPO", "."))
library(data.table)
for (expr in parse(path.expand(Sys.getenv("RIBOCRYPT_RUN_SCRIPT", "~/Desktop/run_ribocrypt.R")))) {
  if (!is.call(expr)) next
  fun <- as.character(expr[[1]])
  if (identical(fun, "RiboCrypt_app")) break
  if (identical(fun, "library") || any(fun == "load_all")) next
  eval(expr)
}
directory <- "/tmp/ribocrypt-account/data"
dir.create(directory, recursive = TRUE, showWarnings = FALSE)
database <- file.path(directory, "access.sqlite")
if (file.exists(database)) unlink(database)
con <- ribocrypt_access_db(database)
exp_dir <- ORFik::config()["exp"]
for (name in c("human_all_merged_l50", "all_samples-Homo_sapiens"))
  ribocrypt_access_dataset(con, name, name, exp_dir, public = TRUE, reference_annotations = TRUE)
base <- read.table(file.path(exp_dir, "human_all_merged_l50.csv"), header = FALSE,
                   sep = ",", fill = TRUE, colClasses = "character", quote = "\"")
for (user in c("alice", "bob")) {
  name <- paste0(user, "_private")
  clone <- base[1:5, ]
  clone[1, 2] <- name
  if (!"Run" %in% clone[4, ]) {
    clone[[ncol(clone) + 1L]] <- ""
    clone[4, ncol(clone)] <- "Run"
  }
  run_column <- which(clone[4, ] == "Run")
  clone[5, run_column] <- paste0("PRIVATE_", toupper(user), "_RUN")
  ORFik::save.experiment(clone, file.path(directory, paste0(name, ".csv")))
  ribocrypt_access_workspace(con, user)
  ribocrypt_access_member(con, Sys.getenv("RIBOCRYPT_TEST_ISSUER", "test-issuer"),
                         Sys.getenv(paste0("RIBOCRYPT_TEST_", toupper(user), "_SUBJECT"), user), user)
  ribocrypt_access_dataset(con, name, name, directory, reference_annotations = TRUE)
  ribocrypt_access_grant(con, user, name, download = user == "bob")
  extra <- copy(m[1])
  extra$Run <- paste0("PRIVATE_", toupper(user), "_RUN")
  extra$sample_title <- paste0("ONLY_VISIBLE_TO_", toupper(user))
  m <- rbind(m, extra, fill = TRUE)
}
DBI::dbDisconnect(con)
expr$access_control <- ribocrypt_access_control(
  database, Sys.getenv("RIBOCRYPT_TEST_ISSUER", "test-issuer"),
  Sys.getenv("RIBOCRYPT_TEST_GATEWAY_SECRET", paste(rep("t", 64), collapse = "")))
expr$init_tab_focus <- "browser"
app <- eval(expr)
shiny::runApp(app, host = "127.0.0.1", port = 7837, launch.browser = FALSE)
