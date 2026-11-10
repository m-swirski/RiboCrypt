#!/usr/bin/env Rscript
# Cache descriptive, default WT/first-condition contrasts. Replication must be
# reviewed in the UI before inference; this script creates no inferential claims.
args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)==1L)
file_arg<-sub('^--file=','',commandArgs()[grepl('^--file=',commandArgs())])
suppressMessages(devtools::load_all(normalizePath(file.path(dirname(file_arg),'../..')),quiet=TRUE))
dir<-normalizePath(args[1]);bundle<-translation_read_bundle(dir)
store<-readRDS(file.path(dir,'counts.rds'));mapping<-as.data.frame(data.table::fread(file.path(dir,'run_map.csv'),colClasses='character'))
p<-translation_pair_store(store,mapping)
if(file.exists(file.path(dir,'pair_metadata.rds'))) {
  qc<-readRDS(file.path(dir,'pair_metadata.rds'))
  p<-translation_attach_pair_qc(p,qc)
  # Default UI uses this lightweight, frozen metadata with the reference cutoff.
  p$units<-translation_filter_pair_qc(qc$units,minimum=.188)
}
cache<-file.path(dir,'fit_cache');dir.create(cache,showWarnings=FALSE)
audit<-list()
for(study in sort(unique(na.omit(p$units$study)))) {
 u0<-p$units[p$units$study==study,,drop=FALSE]
 contexts<-unique(u0$context[!is.na(u0$context)])
 for(ctx in unique(c('__all__',contexts))) {
  u<-if(ctx=='__all__')u0 else u0[!is.na(u0$context)&u0$context==ctx,,drop=FALSE]
  cs<-sort(unique(na.omit(u$condition)));cs<-cs[nzchar(cs)];if(length(cs)<2L)next
  ref<-if('WT'%in%cs)'WT' else cs[1]
  for(cmp in setdiff(cs,ref)) {
   selected<-u[u$condition%in%c(ref,cmp),,drop=FALSE]
   payload<-list(Ribo=p$Ribo[,selected$biosample_id,drop=FALSE],RNA=p$RNA[,selected$biosample_id,drop=FALSE],units=selected,
      reference=ref,comparison=cmp,descriptive=TRUE,cacheable=TRUE,
      reviewed_replication=FALSE,paired=FALSE,batch=NULL,genes=p$genes,count_region=p$count_region)
   payload$source_key<-digest::digest(list(payload,source=bundle$source_run,reference=p$reference,engine=body(translation_fit)),algo='sha256')
   path<-file.path(cache,paste0(payload$source_key,'.rds'));t<-proc.time()[3]
   ans<-tryCatch({if(file.exists(path))readRDS(path) else {z<-translation_fit(payload);saveRDS(z,path);z}},error=function(e)e)
   audit[[length(audit)+1L]]<-data.frame(study=study,context=ctx,reference=ref,comparison=cmp,units=nrow(selected),
       status=if(inherits(ans,'error'))conditionMessage(ans) else 'cached_descriptive',
       features=if(inherits(ans,'error'))NA_integer_ else nrow(ans$table),seconds=proc.time()[3]-t,cache_key=payload$source_key)
  }
 }
}
data.table::fwrite(data.table::rbindlist(audit),file.path(dir,'descriptive_precompute_audit.csv'))
cat('Cached',sum(vapply(audit,function(x)x$status=='cached_descriptive',logical(1))),'descriptive contrasts.\n')
