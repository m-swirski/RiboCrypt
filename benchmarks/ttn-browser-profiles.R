# Run from the checkout: OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript benchmarks/ttn-browser-profiles.R
# Uses local public all_merged human data, not the live application or a browser.
suppressPackageStartupMessages({library(ORFik);library(data.table);library(BiocParallel)})
data.table::setDTthreads(1)
devtools::load_all('.',quiet=TRUE)
df <- ORFik::read.experiment('all_merged-Homo_sapiens',validate=FALSE)
grl <- ORFik::loadRegion(df,'mrna')['ENST00000589042']
reads <- RiboCrypt:::get_track_paths(df)[[1]]
stopifnot(length(grl)==1L, length(grl[[1]])==363L)
old <- function(n) lapply(seq_len(n),function(i) RiboCrypt:::getProfileWrapper(grl,reads,TRUE,9,FALSE,type='area'))
new <- function(n) RiboCrypt:::build_browser_profiles(grl,rep(list(reads),n),9,'mean','area','all',rep(TRUE,n),rep(FALSE,n),BiocParallel::SerialParam())
measure <- function(fun,n,reps=4L){fun(n);vapply(seq_len(reps),function(i){start<-proc.time()[['elapsed']];invisible(fun(n));proc.time()[['elapsed']]-start},numeric(1))}
for (n in c(1L,4L)) {
 a<-old(n);b<-new(n);stopifnot(identical(a,b));rm(a,b)
 for (round in 1:2) {a<-measure(old,n,2);b<-measure(new,n,2);cat('tracks',n,'round',round,'baseline',paste(a,collapse=','),'optimized',paste(b,collapse=','),'\n')}
}
cat('Version R',R.version.string,'data.table',as.character(packageVersion('data.table')),'ORFik',as.character(packageVersion('ORFik')),'\n')
