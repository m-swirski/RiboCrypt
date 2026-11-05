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
for(i in 1:3){cat('TTN BACKEND TRIAL',i,'\n');start<-proc.time()[['elapsed']];bottom<-RiboCrypt:::bottom_panel_shiny(fn,templates=templates);p<-RiboCrypt:::browser_track_panel_shiny(fn,bottom,session,templates=templates);cat('BOTTOM+TRACK+MERGE',proc.time()[['elapsed']]-start,'s\n')}

# Full widget JSON preparation; includes plotly_build, not web transfer or browser draw.
for(i in 1:4){start<-proc.time()[['elapsed']];z<-plotly::plotly_json(p,jsonedit=FALSE,pretty=FALSE);cat('FULL PLOTLY JSON',i,proc.time()[['elapsed']]-start,'s, bytes',nchar(z,type='bytes'),'\n')}
Rprof('/tmp/ribocrypt-ttn-backend-Rprof.out',interval=.01,memory.profiling=TRUE,gc.profiling=TRUE)
for(i in 1:6){bottom<-RiboCrypt:::bottom_panel_shiny(fn,templates=templates);p<-RiboCrypt:::browser_track_panel_shiny(fn,bottom,session,templates=templates);z<-plotly::plotly_json(p,jsonedit=FALSE,pretty=FALSE)}
Rprof(NULL);cat('RPROF TOP TOTAL\n');print(head(summaryRprof('/tmp/ribocrypt-ttn-backend-Rprof.out',memory='both')$by.total,30));cat('RPROF TOP SELF\n');print(head(summaryRprof('/tmp/ribocrypt-ttn-backend-Rprof.out',memory='both')$by.self,20))
