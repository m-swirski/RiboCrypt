#!/usr/bin/env Rscript
# Read a completed matched TE run and export a new, opt-in display release.
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=2L)stop('Usage: Rscript prepare_translation_atlas.R completed_TE_run new_output_directory')
run<-normalizePath(args[1]);out<-path.expand(args[2]);if(dir.exists(out))stop('Preserve existing releases; choose a new output directory.')
suppressPackageStartupMessages({library(data.table);library(digest)})
valid<-fread(file.path(run,'run_validation.csv'));stopifnot(all(valid$passed))
dir.create(out,recursive=TRUE)
p<-readRDS(file.path(run,'fig4/functional_profiles.rds'))
native<-readRDS(file.path(run,'fig4/native_canonical_counts.rds'))
libs<-native$libraries;pm<-p$pairs
R<-native$Ribo[,pm$ribo_library_key,drop=FALSE];A<-native$RNA[,pm$rna_library_key,drop=FALSE]
stopifnot(identical(rownames(R),rownames(A)),!anyDuplicated(colnames(R)),!anyDuplicated(colnames(A)))
ids<-pm$pair_id
units<-data.frame(biosample_id=ids,ribo_library_id=pm$ribo_library_key,rna_library_id=pm$rna_library_key,
                 condition=pm$ribo_CONDITION,study=pm$ribo_BioProject,replicate=pm$ribo_REPLICATE,
                 match_level=pm$biological_specimen_pair_status,context=pm$profile_context,
                 ribo_experiment=paste0(pm$ribo_BioProject,'-homo_sapiens'),
                 rna_experiment=paste0(pm$rna_BioProject,'-homo_sapiens_RNA-seq'))
# Resolve true experiment paths from the verified native Run inventory, never
# assume that BioProject and processing study accession are interchangeable.
inv<-fread(file.path(run,'fig4/native_run_inventory.csv'))
members<-rbindlist(lapply(seq_len(nrow(units)),function(i) {
 rbindlist(lapply(c('RIBO','RNA'),function(assay) {
  key<-if(assay=='RIBO')pm$ribo_library_key[i] else pm$rna_library_key[i]
  lib<-libs[library_key==key];stopifnot(nrow(lib)==1L)
  rr<-strsplit(lib$run_ids,';',fixed=TRUE)[[1]]
  data.table(Run=rr,library_id=key,assay=assay,biosample_id=ids[i])
 }))
}))
stopifnot(!anyDuplicated(members$Run))
for(i in seq_len(nrow(units)))for(assay in c('RIBO','RNA')) {
 key<-if(assay=='RIBO')units$ribo_library_id[i] else units$rna_library_id[i]
 paths<-unique(inv[library_key==key & inv$assay==assay]$experiment_file)
 stopifnot(length(paths)==1L)
 units[[if(assay=='RIBO')'ribo_experiment' else 'rna_experiment']][i]<-sub('\\.csv$','',basename(paths))
}
# Author-supported matching remains a proposal; it is not upgraded to same material.
map<-merge(members[,.(Run,biosample_id)],as.data.table(units)[,.(biosample_id,condition,study,replicate,match_level)],by='biosample_id',sort=FALSE)
reference<-digest(file.path(run,'fig4/gene_transcript_reference.csv'),file=TRUE,algo='sha256')
store<-list(version=1L,Ribo=R,RNA=A,run_members=as.data.frame(members[,.(Run,library_id,assay)]),units=units,genes=as.data.frame(native$genes),count_region='Ribo CDS / RNA CDS, canonical transcript per gene',reference=reference,source_run=run)
saveRDS(store,file.path(out,'counts.rds'))
fwrite(map,file.path(out,'run_map.csv'))
extras<-list(terms=fread(file.path(run,'fig4/GO_functional_AUROC.csv')),
             candidates=fread(file.path(run,'fig5/current_TE_specific_shortlist.csv')),
             retained=fread(file.path(run,'model_screen/all_models_retained_candidates.csv')),
             memberships=fread(file.path(run,'fig4/measured_gene_GO_memberships.csv')))
experiments<-unique(c(units$ribo_experiment,units$rna_experiment))
bundle<-c(list(version=1L,species='Homo sapiens',TE=p$profiles$TE,RNA=p$profiles$RNA,genes=as.data.frame(p$genes),contexts=as.data.frame(p$contexts),source_run=run,source_experiments=experiments,visibility='public',pooling='Mean within study/context, then equal mean across studies; contexts may mix treatments',TE_definition='Versioned Cenik-style ILR regression / CLR reconstructed residuals; gene order and universe frozen',reference=reference),extras)
saveRDS(bundle,file.path(out,'atlas.rds'))
# The full all-gene matrix is deliberately not serialized for browsing.
x<-bundle$TE;xc<-x-rowMeans(x);v<-rowSums(xc^2)
nei<-rbindlist(lapply(seq_len(nrow(x)),function(i) {
 rho<-2*drop(xc%*%xc[i,])/(v+v[i]);rho[i]<-NA_real_;rho[!is.finite(rho)]<-NA_real_
 ix<-head(order(abs(rho),decreasing=TRUE,na.last=NA),50)
 data.table(query=rownames(x)[i],gene_id=rownames(x)[ix],TE_rho=rho[ix])
}))
saveRDS(nei,file.path(out,'top_neighbours.rds'))
sources<-c('run_validation.csv','fig4/functional_profiles.rds','fig4/native_canonical_counts.rds','fig4/GO_functional_AUROC.csv','fig4/measured_gene_GO_memberships.csv','fig5/current_TE_specific_shortlist.csv','model_screen/all_models_retained_candidates.csv')
manifest<-data.table(source=file.path(run,sources),sha256=vapply(file.path(run,sources),digest,character(1),file=TRUE,algo='sha256'))
fwrite(manifest,file.path(out,'source_manifest.csv'))
files<-list.files(out,full.names=TRUE);fwrite(data.table(file=basename(files),bytes=file.info(files)$size,sha256=vapply(files,digest,character(1),file=TRUE,algo='sha256')),file.path(out,'output_manifest.csv'))
writeLines(c('Opt-in public Translation release. Curated mapping has not established same material.',paste('Source:',run),paste('Units:',nrow(units)),paste('Genes:',nrow(x)),paste('Contexts:',ncol(x))),file.path(out,'README.txt'))
cat('Prepared',out,'\n')
