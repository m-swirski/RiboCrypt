# Compare current vs clean-HEAD sequence payload code in one R process.
# git show 3c72c86:R/js_helpers.R > /tmp/ribocrypt_old_js_helpers.R
suppressPackageStartupMessages({library(ORFik);library(data.table)});setDTthreads(1);devtools::load_all('.',quiet=TRUE)
old<-new.env(parent=asNamespace('RiboCrypt'));sys.source('/tmp/ribocrypt_old_js_helpers.R',envir=old)
df<-read.experiment('all_merged-Homo_sapiens',validate=FALSE);gr<-loadRegion(df,'mrna')['ENST00000589042'];s<-RiboCrypt:::browser_tx_seqs_getSeq(findFa(df),gr)
f<-list(old=function()old$fetch_JS_seq(s,1,250,'one_letter','selectedRegion','R'),new=function()RiboCrypt:::fetch_JS_seq(s,1,250,'one_letter','selectedRegion','R'))
for(n in names(f)){x<-f[[n]]();cat('PAYLOAD',n,nchar(jsonlite::toJSON(x,auto_unbox=TRUE,null='null'),type='bytes'),'bytes\n')}
for(i in 1:8)for(n in if(i%%2)c('old','new')else c('new','old')){t<-proc.time()[['elapsed']];f[[n]]();cat('TIME',i,n,proc.time()[['elapsed']]-t,'\n')}
