# Run after devtools::load_all() on the home-selector branch.
# These are local-server paths and can be adjusted for another installation.
library(data.table)
all_exp <- ORFik::list.experiments(validate = FALSE)
# Keep per-modality collections for the Browser and Home. MegaBrowser's
# collection list remains restricted to the existing canonical collections.
all_exp_meta <- all_exp[grepl("^all_samples-[A-Za-z]+_[a-z]+$", name)]
keep_samples <- grepl("^all_samples-[A-Za-z]+_[a-z]+(_(RNA|SSU|disome|CAGE))?$", all_exp$name, ignore.case=TRUE)
all_exp <- all_exp[!grepl("all_samples", name) | keep_samples]
metadata <- fread("~/livemount/Bio_data/NGS_pipeline/metadata_done_samples_extended.csv")
app <- RiboCrypt::RiboCrypt_app(
  validate.experiments = FALSE, all_exp = all_exp, all_exp_meta = all_exp_meta,
  metadata = metadata,
  browser_options = c(default_experiment = "all_merged-Homo_sapiens",
                      default_experiment_meta = "all_samples-Homo_sapiens"),
  home_html = system.file("home", "home.html", package = "RiboCrypt"),
  init_tab_focus = "home"
)
shiny::runApp(app)
