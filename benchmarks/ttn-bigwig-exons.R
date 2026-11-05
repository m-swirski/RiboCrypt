suppressPackageStartupMessages({library(ORFik);library(data.table);library(GenomicRanges);library(rtracklayer)})
setDTthreads(1);df<-read.experiment('all_merged-Homo_sapiens',validate=FALSE);gr<-loadRegion(df,'mrna')['ENST00000589042'];read<-ORFik::filepath(df,'bigwig',suffix_stem=c('_pshifted',''))[[1]]
width0<-sum(width(unlist(gr))); start0<-178525989L; long<-GRangesList(tx=GRanges('2',IRanges(start0,width=width0),strand='-')); splitlen<-rep(width0 %/% 363L,363);splitlen[1:(width0%%363L)]<-splitlen[1:(width0%%363L)]+1L; starts<-start0+c(0L,head(cumsum(splitlen),-1L));split<-GRangesList(tx=GRanges('2',IRanges(starts,width=splitlen),strand='-'));split<-endoapply(split,rev)
bench<-function(expr,n=6){force(expr); out<-numeric(n);for(i in seq_len(n)){st<-proc.time()[['elapsed']];force(expr());out[i]<-proc.time()[['elapsed']]-st};out}
for(name in c('ttn','contiguous','synthetic363')){x<-switch(name,ttn=gr,contiguous=long,synthetic363=split);fun<-function()ORFik::coveragePerTiling(x,read,as.data.table=TRUE,withFrames=TRUE,is.sorted=TRUE);fun();cat(name,'exons',length(x[[1]]),'width',sum(width(unlist(x))),'seconds',paste(bench(fun),collapse=','),'\n')}
rl<-ranges(gr);names(rl)<-as.character(seqnames(gr[[1]])[1]);rl<-rl[1];temp<-rtracklayer::import.bw(read[[2]],as='NumericList',which=rl)
old<-function(){v<-rev(temp);unlist(relist(rev(unlist(v)),skeleton=temp),use.names=FALSE)}
new<-function(){v<-unlist(temp,use.names=FALSE);n<-lengths(temp);v[sequence(n,from=cumsum(n),by=-1L)]}
stopifnot(identical(old(),new()));cat('reverse old',paste(bench(old,10),collapse=','),'new',paste(bench(new,10),collapse=','),'\n')
