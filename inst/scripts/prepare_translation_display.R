#!/usr/bin/env Rscript
# Add small display/QC sidecars to a new immutable release; keep source release intact.
args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)==2L)
source<-normalizePath(args[1]);out<-path.expand(args[2]);if(dir.exists(out))stop('Choose a new display release directory.')
suppressPackageStartupMessages({library(data.table);library(digest)})
dir.create(out,recursive=TRUE)
# Byte-identical read-only resources can share disk blocks on the data mount.
for(name in c('atlas.rds','counts.rds','run_map.csv','top_neighbours.rds','source_manifest.csv')) {
 if(!file.link(file.path(source,name),file.path(out,name)))stopifnot(file.copy(file.path(source,name),file.path(out,name)))
}
b<-readRDS(file.path(source,'atlas.rds'));s<-readRDS(file.path(source,'counts.rds'));run<-b$source_run
qa<-as.data.frame(fread(file.path(run,'fig4/pair_QC_and_cohort_audit.csv')))
units<-s$units;ix<-match(units$biosample_id,qa$pair_id);stopifnot(!anyNA(ix))
units$matched_R2<-qa$matched_R2[ix];units$matching_score<-qa$matching_score[ix];units$Ribo_CDS_fraction_proxy<-qa$ribo_CDS_fraction_proxy[ix]
units$qc_ribo_library_id<-qa$ribo_library_key[ix];units$qc_rna_library_id<-qa$rna_library_key[ix]
frozen<-as.data.frame(fread(file.path(run,'fig4/gene_filter_audit.csv')))
qc<-list(units=units,run_members=s$run_members,reference=s$reference,
   feature_ids=frozen$gene_id[frozen$expression_filter],other_ids=frozen$gene_id[!frozen$replication_dependent_histone_proxy & !frozen$expression_filter],
   definition='Squared Pearson correlation across matched RNA/Ribo CLR-count components; zeros floored to one; frozen expressed genes plus others component',reference_cutoff=.188)
saveRDS(qc,file.path(out,'pair_metadata.rds'))
core<-b[setdiff(names(b),c('terms','candidates','retained','memberships'))]
saveRDS(core,file.path(out,'atlas_core.rds'),compress=FALSE)
go<-b[c('terms','candidates','retained','memberships')]
# The display needs one gene/term membership, not repeated annotation evidence rows.
keep<-go$memberships$GO_ID %in% go$terms$GO_ID
go$memberships<-unique(go$memberships[keep,c('GO_ID','gene_id','ontology','GO_name'),with=FALSE])
go$memberships$symbol<-b$genes$symbol[match(go$memberships$gene_id,b$genes$gene_id)]
con<-gzfile(file.path(out,'atlas_go.rds'),'wb',compression=1);saveRDS(go,con);close(con)
metrics<-fread(file.path(run,'fig4/cohort_summary.csv'))
metric<-function(name) metrics$value[match(name,metrics$metric)]
# Minimal metadata, never raw counts, backs the landing page.
q<-qa[,c('pair_id','ribo_BioProject','profile_context','both_complete_counts','matched_R2','matching_score','ribo_CDS_fraction_proxy','included_in_primary','exclusion_reason')]
cover<-as.data.frame(b$contexts)
by_study<-as.data.frame(as.data.table(units)[,.(pairs=.N),by=study][order(-pairs)])
exclusions<-as.data.frame(as.data.table(q)[,.(pairs=.N),by=exclusion_reason][order(-pairs)])
overview<-list(species=b$species,counts=c(curated_pairs=metric('strict_metadata_pairs'),count_ready_pairs=metric('pairs_complete_native_counts'),R2_passing_pairs=metric('pairs_R2_reference_pass'),CDS_and_R2_passing_pairs=metric('pairs_R2_and_CDS_proxy_pass'),atlas_pairs=metric('primary_pairs'),contexts=metric('primary_contexts'),studies=metric('primary_studies'),genes=metric('functional_genes')),
  pairs=q,contexts=cover,studies=by_study,exclusions=exclusions,reference_cutoff=.188,source_run=run,
  source_experiments=unique(as.character(fread(file.path(run,'discovery/eligible_experiments.csv'))$experiment)))
saveRDS(overview,file.path(out,'overview.rds'),compress=FALSE)
header<-list(version=1L,species=b$species,genes=nrow(b$TE),contexts=ncol(b$TE),pairs=nrow(units),source_run=run,source_experiments=b$source_experiments)
saveRDS(header,file.path(out,'dataset.rds'),compress=FALSE)
jsonlite::write_json(list(source_release=source,source_run=run,metric_definition=qc$definition,reference_cutoff=.188,prepared=format(Sys.time(),tz='UTC')),file.path(out,'display_provenance.json'),auto_unbox=TRUE,pretty=TRUE)
files<-list.files(out,full.names=TRUE);fwrite(data.table(file=basename(files),bytes=file.info(files)$size,sha256=vapply(files,digest,character(1),file=TRUE,algo='sha256')),file.path(out,'output_manifest.csv'))
cat('Prepared display release:',out,'\n')
