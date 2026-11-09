# Run after loading/installing the updated RiboCrypt package. The script only
# builds default 2D tissue/cell-line caches; missing/stale inputs are reported.
suppressPackageStartupMessages({library(data.table);library(ORFik);library(RiboCrypt)})
exps <- ORFik::list.experiments(validate=FALSE,pattern='all_samples',libtypeExclusive='RFP')
exps <- exps[grepl('^all_samples-[A-Za-z]+_[a-z]+$',name) & name != 'all_samples-all_species']
output <- Sys.getenv('RIBOCRYPT_UMAP_REPORT','observatory_umap_cache_status.csv')
status <- list()
for (id in exps$name) {
  t <- Sys.time()
  result <- tryCatch({
    e <- ORFik::read.experiment(id,validate=FALSE)
    source <- file.path(ORFik::refFolder(e),'UMAP','UMAP_by_gene_counts.fst')
    if (!file.exists(source)) 'no_umap' else {RiboCrypt::precompute_observatory_umap(e); 'done'}
  },error=function(e) paste0('ERROR: ',conditionMessage(e)))
  status[[length(status)+1L]] <- data.table(experiment=id,status=result,seconds=as.numeric(difftime(Sys.time(),t,units='secs')))
  fwrite(rbindlist(status),output)
  message(id,': ',result)
}
