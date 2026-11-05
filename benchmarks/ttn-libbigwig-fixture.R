# Research prototype regression fixtures; not a production acceptance suite.
suppressPackageStartupMessages({library(ORFik);library(data.table);library(GenomicRanges);library(rtracklayer);library(Rcpp)})
setDTthreads(1);prefix<-Sys.getenv('LIBBIGWIG_PREFIX');stopifnot(nzchar(prefix));Sys.setenv(PKG_CPPFLAGS=paste0('-I',prefix,'/include'),PKG_LIBS=paste0('-L',prefix,'/lib -lBigWig'))
Rcpp::sourceCpp('benchmarks/ttn-libbigwig-prototype.cpp',rebuild=TRUE)
bw<-tempfile(fileext='.bw');on.exit(unlink(bw)); chr<-Seqinfo::Seqinfo('1',seqlengths=50000)
for (scores in list(c(0,3,4,7,1),c(0,1.5,-2.5,0.125,10))) {
 signal<-GRanges('1',IRanges(c(100,500,1000,3000,10000),width=c(30,100,120,50,20)),score=scores,seqinfo=chr);export(signal,bw,format='bigWig')
 for(minus in c(FALSE,TRUE)) {
  ex<-GRanges('1',IRanges(start=c(80,100,115,490,510,1100,2000,2999,9995,12000,12000),width=c(100,20,1,80,90,50,30,51,40,50,50)),strand=if(minus)'-' else '+',seqinfo=chr)
  if(minus)ex<-rev(ex);gr<-GRangesList(tx=ex);base<-coveragePerTiling(gr,c(bw,bw),as.data.table=TRUE,withFrames=TRUE,is.sorted=TRUE)$count
  for(g in c(-1L,0L,10L,500L,50000L)){val<-libbigwig_exons(bw,'1',start(ex),end(ex),minus,g);if(!identical(as.numeric(base),as.numeric(val)))stop('fixture mismatch');cat('PASS scores=',paste(scores,collapse=','),'minus=',minus,'gap=',g,'\n')}
 }
}
