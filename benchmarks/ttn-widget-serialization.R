# Single-track real TTN backend (not browser/network latency).
suppressPackageStartupMessages({library(ORFik);library(data.table)})
setDTthreads(1);devtools::load_all('.',quiet=TRUE);options(RiboCrypt.browser_profile=TRUE)
df<-read.experiment('all_merged-Homo_sapiens',validate=FALSE);tx<-loadRegion(df,'mrna');cds<-loadRegion(df,'cds');id<-'ENST00000589042'
display<-RiboCrypt:::observed_tx_annotation(id,function()tx)
a<-RiboCrypt:::observed_cds_annotation_internal(id,cds,FALSE)
ta<-RiboCrypt:::observed_cds_annotation_internal(id,tx,FALSE)
controls<-list(dff=df,display_region=display,annotation=a,extendLeaders=0,extendTrailers=0,viewMode=FALSE,custom_sequence='',customRegions=GenomicRanges::GRangesList(),tx_annotation=ta,collapsed_introns_width=0,frame_colors='R',phyloP=FALSE,mapability=FALSE,is_cellphone=FALSE,reads=RiboCrypt:::get_track_paths(df),withFrames=TRUE,frames_type='area',colors=NULL,kmerLength=9,log_scale=FALSE,summary_track=FALSE,summary_track_type='area',export_format='svg',zoom_range=numeric(),frames_subset='all',y_range=NULL)
fn<-function()controls;session<-list(ns=identity)
templates<-list(nt_seq_panel_plotly=RiboCrypt:::ntSeqPanelPlotlyTemplate(),aa_seq_panel_plotly=RiboCrypt:::aaSeqPanelPlotlyTemplate(),gene_model_panel_plotly=RiboCrypt:::geneModelPanelPlotlyTemplate(),cov_panel_area_plotly=RiboCrypt:::covPanelAreaPlotlyTemplate())
bottom<-RiboCrypt:::bottom_panel_shiny(fn,templates=templates)
p<-RiboCrypt:::browser_track_panel_shiny(fn,bottom,session,templates=templates)
oldenv<-new.env(parent=asNamespace('RiboCrypt'));sys.source('/tmp/ribocrypt_old_js_helpers.R',oldenv)
oldpayload<-oldenv$fetch_JS_seq(bottom$target_seq,0,250,'one_letter','selectedRegion','R')
newpayload<-RiboCrypt:::fetch_JS_seq(bottom$target_seq,0,250,'one_letter','selectedRegion','R')
# Reuse the real data and hooks, replacing only the sequence hook payload.
q<-plotly::plotly_build(p); old<-new<-q
old$jsHooks$render[[1]]$data<-oldpayload;new$jsHooks$render[[1]]$data<-newpayload
obj<-list(old=list(x=old$x,evals=NULL,jsHooks=old$jsHooks),new=list(x=new$x,evals=NULL,jsHooks=new$jsHooks))
for(n in names(obj)){z<-htmlwidgets:::toJSON(obj[[n]]);cat('WIDGET BYTES',n,nchar(z,type='bytes'),'\n')}
for(i in 1:8)for(n in if(i%%2)c('old','new') else c('new','old')){t<-proc.time()[['elapsed']];z<-htmlwidgets:::toJSON(obj[[n]]);cat('SERIALIZE',i,n,proc.time()[['elapsed']]-t,'\n')}
