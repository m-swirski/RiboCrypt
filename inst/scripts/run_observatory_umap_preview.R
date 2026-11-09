# Run after devtools::load_all() on this worktree. Local server paths only.
library(data.table)
all_exp <- ORFik::list.experiments(validate=FALSE)
all_exp_meta <- all_exp[grepl('^all_samples-[A-Za-z]+_[a-z]+$',name)]
all_exp <- rbindlist(list(all_exp[!grepl('all_samples',name)],all_exp_meta))
metadata <- fread('~/livemount/Bio_data/NGS_pipeline/metadata_done_samples_extended.csv')
app <- RiboCrypt::RiboCrypt_app(validate.experiments=FALSE,all_exp=all_exp,all_exp_meta=all_exp_meta,
  metadata=metadata,init_tab_focus='Observatory',
  browser_options=c(default_experiment='all_merged-Homo_sapiens',
                    default_experiment_meta='all_samples-Homo_sapiens'))
shiny::runApp(app)
