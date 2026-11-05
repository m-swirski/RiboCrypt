# Reproduce frame preparation benchmark after devtools::load_all(".").
suppressPackageStartupMessages(library(data.table));setDTthreads(1);devtools::load_all('.',quiet=TRUE)
set.seed(4);for(size in list(c(1434L,1000L),c(2000L,1500L))){
 mat<-matrix(rexp(prod(size)),nrow=size[1]);measure<-function(f){f();vapply(1:5,function(i){t<-proc.time()[['elapsed']];f();proc.time()[['elapsed']]-t},numeric(1))}
 original<-function(){finite<-is.finite(mat);limits<-range(mat[finite]);if(limits[1]>0)limits[1]<-0;x<-pmin(pmax((mat-limits[1])/diff(limits),0),1)^0.5*(1-1e-7);z<-sweep(x,2,(seq_len(ncol(mat))-1L)%%3L,'+');z[!finite]<-NA_real_;z}
 fast<-function()RiboCrypt:::megabrowser_frame_encode(mat)
 stopifnot(isTRUE(all.equal(original(),fast(),tolerance=0)));cat('size',paste(size,collapse='x'),'original',paste(measure(original),collapse=','),'optimized',paste(measure(fast),collapse=','),'\n')
}
