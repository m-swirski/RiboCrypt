# Validate a prepared release, exact pair/QC identities, caches and lazy plots.
args<-commandArgs(TRUE);directory<-normalizePath(args[1]);evidence<-normalizePath(args[2])
file_arg<-sub('^--file=','',commandArgs()[grepl('^--file=',commandArgs())])
suppressMessages(devtools::load_all(normalizePath(file.path(dirname(file_arg),'../..')),quiet=TRUE))
b<-RiboCrypt:::translation_read_bundle(directory);store<-readRDS(file.path(directory,'counts.rds'));mapping<-as.data.frame(data.table::fread(file.path(directory,'run_map.csv'),colClasses='character'))
p<-RiboCrypt:::translation_pair_store(store,mapping);qc<-readRDS(file.path(directory,'pair_metadata.rds'));p<-RiboCrypt:::translation_attach_pair_qc(p,qc)
stopifnot(identical(dimnames(p$Ribo),dimnames(p$RNA)),all(is.finite(p$Ribo)),all(is.finite(p$RNA)),all(p$Ribo>=0),all(p$RNA>=0),!anyDuplicated(p$units$biosample_id),setequal(p$units$biosample_id,store$units$biosample_id))
o<-readRDS(file.path(directory,'overview.rds'));stopifnot(o$counts[['atlas_pairs']]==nrow(p$units),o$counts[['genes']]==nrow(b$TE),o$counts[['contexts']]==ncol(b$TE),all(p$units$biosample_id %in% o$pairs$pair_id))
pre<-readRDS(file.path(directory,'top_neighbours.rds'));genes<-unique(c(b$genes$gene_id[match(c('ATF4','EIF5A','DOHH'),b$genes$symbol)],head(b$genes$gene_id,2)));genes<-genes[!is.na(genes)]
for(g in genes) {
 z<-RiboCrypt:::translation_neighbours(b,g,n=30,precomputed=pre);raw<-RiboCrypt:::translation_neighbours(b,g,n=30)
 stopifnot(identical(z$gene_id,raw$gene_id),isTRUE(all.equal(z$TE_rho,raw$TE_rho,tolerance=1e-12)),isTRUE(all.equal(z$RNA_rho,raw$RNA_rho,tolerance=1e-12)))
}
audit<-data.table::fread(file.path(directory,'descriptive_precompute_audit.csv'));stopifnot(nrow(audit)>0,all(audit$status=='cached_descriptive'),all(file.exists(file.path(directory,'fit_cache',paste0(audit$cache_key,'.rds')))))
rm(list=ls(RiboCrypt:::.translation_cache),envir=RiboCrypt:::.translation_cache)
times<-list();g<-genes[1]
shiny::testServer(RiboCrypt:::translation_server,args=list(browser_options=c(translation_data_dir=directory)),{
 t<-proc.time()[3];session$setInputs(species='Homo sapiens',mode='Species overview');session$flushReact()
 stopifnot(grepl('Curated assay pairs',output$overview_summary$html,fixed=TRUE),nzchar(output$overview_qc),nzchar(output$overview_coverage))
 stopifnot(!any(grepl('counts.rds|atlas.rds|atlas_core.rds|top_neighbours.rds',ls(RiboCrypt:::.translation_cache))))
 times$overview_seconds<<-unname(proc.time()[3]-t)
 t<-proc.time()[3];session$setInputs(mode='Which genes move together?',gene=g,contexts=colnames(b$TE),n_neighbours=30);session$flushReact()
 stopifnot(nzchar(output$profiles),nzchar(output$network));times$gene_plots_seconds<<-unname(proc.time()[3]-t)
 t<-proc.time()[3];session$setInputs(mode='What functions are suggested?',pathway=b$terms$GO_ID[1],broad=FALSE);session$flushReact()
 stopifnot(nzchar(output$term_performance));times$GO_plot_seconds<<-unname(proc.time()[3]-t)
})
jsonlite::write_json(list(directory=directory,source_run=b$source_run,pairs=nrow(p$units),genes=nrow(b$TE),contexts=ncol(b$TE),cached_contrasts=nrow(audit),plots=times,validation='passed'),file.path(evidence,'atlas_validation.json'),auto_unbox=TRUE,pretty=TRUE)
cat('PASS: real pairing, QC identities, cached signed neighbours, contrast cache and lazy overview/gene/GO rendering\n')
