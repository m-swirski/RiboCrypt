# Single-track real TTN backend (not browser/network latency).
suppressPackageStartupMessages({library(ORFik);library(data.table)})
setDTthreads(1);devtools::load_all('.',quiet=TRUE);options(RiboCrypt.browser_profile=TRUE)
df<-read.experiment('all_merged-Homo_sapiens',validate=FALSE);tx<-loadRegion(df,'mrna');cds<-loadRegion(df,'cds');id<-'ENST00000589042'
display<-RiboCrypt:::observed_tx_annotation(id,function()tx)
a<-RiboCrypt:::observed_cds_annotation_internal(id,cds,FALSE)
ta<-RiboCrypt:::observed_cds_annotation_internal(id,tx,FALSE)
controls<-list(dff=df,display_region=display,annotation=a,extendLeaders=0,extendTrailers=0,viewMode=FALSE,custom_sequence='',customRegions=GenomicRanges::GRangesList(),tx_annotation=ta,collapsed_introns_width=0,frame_colors='R',phyloP=FALSE,mapability=FALSE,is_cellphone=FALSE,reads=RiboCrypt:::get_track_paths(df),withFrames=TRUE,frames_type='area',colors=NULL,kmerLength=9,log_scale=FALSE,summary_track=FALSE,summary_track_type='area',export_format='svg',zoom_range=numeric(),frames_subset='all',y_range=NULL)
fn<-function()controls;session<-list(ns=identity)
templates<-list(nt_seq_panel_plotly=RiboCrypt:::ntSeqPanelPlotlyTemplate(),aa_seq_panel_plotly=RiboCrypt:::aaSeqPanelPlotlyTemplate(),gene_model_panel_plotly=RiboCrypt:::geneModelPanelPlotlyTemplate(),cov_panel_area_plotly=RiboCrypt:::covPanelAreaPlotlyTemplate(),cov_panel_with_frames_plotly=RiboCrypt:::covPanelWithFramesPlotlyTemplate(),cov_panel_columns_plotly=RiboCrypt:::covPanelColumnsPlotlyTemplate())
for(mode in c('lines','columns','area')){controls$frames_type<-mode; controls$summary_track_type<-mode;cat('MODE',mode,'\n');for(i in 1:2){t<-proc.time()[['elapsed']];bottom<-RiboCrypt:::bottom_panel_shiny(fn,templates=templates);p<-RiboCrypt:::browser_track_panel_shiny(fn,bottom,session,templates=templates);cat('BACKEND',mode,i,proc.time()[['elapsed']]-t,'\n');};htmlwidgets::saveWidget(p,paste0('/tmp/ribocrypt_ttn_',mode,'.html'),selfcontained=FALSE);z<-plotly::plotly_json(p,jsonedit=FALSE,pretty=FALSE);cat('JSON',mode,nchar(z,type='bytes'),'\n')}
