translation_ui <- function(id) {
  ns<-shiny::NS(id)
  shiny::tabPanel('Translation',icon=shiny::icon('layer-group'),
    shiny::tags$style('.translation-unit-picker .selectize-input, .translation-context-picker .selectize-input {max-height:160px; overflow-y:auto;} .translation-unit-picker .selectize-dropdown-content {max-height:220px;} .translation-overview-cards {display:flex;flex-wrap:wrap;gap:12px;margin:12px 0;} .translation-overview-card {padding:12px 18px;border:1px solid #ddd;border-radius:6px;min-width:150px;} .translation-overview-card strong {display:block;font-size:25px;}'),
    shiny::sidebarLayout(shiny::sidebarPanel(
      shiny::selectInput(ns('species'),'Species',choices=NULL),
      shiny::uiOutput(ns('status')),
      shiny::conditionalPanel(condition="input.mode !== 'Species overview'",ns=ns,shiny::selectizeInput(ns('gene'),'Gene',choices=NULL)),
      shiny::conditionalPanel(condition="input.mode === 'What changed?'",ns=ns,
      shiny::tags$hr(),
      shiny::selectInput(ns('study'),'Study for differential analysis',choices=NULL),
      shiny::selectInput(ns('cohort_context'),'Cell/tissue context for differential analysis',choices=NULL),
      shiny::checkboxInput(ns('pair_qc'),'Filter pairs by RNA–Ribo correlation QC',TRUE),
      shiny::conditionalPanel(condition='input.pair_qc',ns=ns,
        shiny::numericInput(ns('pair_r2_min'),'Minimum matched RNA–Ribo R² (CLR counts)',value=.188,min=0,max=1,step=.01),
        shiny::checkboxInput(ns('include_unscored'),'Include pairs without a QC score',FALSE)),
      shiny::checkboxInput(ns('all_units'),'Use all eligible matched biosamples',TRUE),
      shiny::uiOutput(ns('unit_count')),
      shiny::conditionalPanel(condition='!input.all_units',ns=ns,
        shiny::div(class='translation-unit-picker',shiny::selectizeInput(ns('units'),'Matched biosamples (keep biological replicates)',choices=NULL,multiple=TRUE,options=list(plugins=list('remove_button'),maxOptions=100)))),
      shiny::selectInput(ns('reference'),'Reference (Group 1)',choices=NULL),
      shiny::selectInput(ns('comparison'),'Comparison (Group 2)',choices=NULL),
      shiny::fileInput(ns('mapping'),'Reviewed Run → biosample mapping (CSV)',accept='.csv'),
      shiny::checkboxInput(ns('descriptive'),'Descriptive effects only',TRUE),
      shiny::checkboxInput(ns('replication'),'Independent biological replicates reviewed',FALSE),
      shiny::checkboxInput(ns('paired'),'RNA/Ribo originate from the same biological material',FALSE),
      shiny::actionButton(ns('fit'),'Run differential TE'),
      shiny::tags$hr(),
      shiny::sliderInput(ns('fdr'),'Adjusted FDR',min=0,max=0.2,value=.05,step=.005),
      shiny::numericInput(ns('effect'),'Minimum |ΔTE| (log₂)',value=.5,min=0,step=.1),
      shiny::downloadButton(ns('download'),'Download differential results')),
      shiny::conditionalPanel(condition="input.mode !== 'Species overview'",ns=ns,shiny::tags$details(shiny::tags$summary('Atlas settings'),
        shiny::div(class='translation-context-picker',shiny::selectizeInput(ns('contexts'),'Atlas contexts',choices=NULL,multiple=TRUE)),
        shiny::numericInput(ns('n_neighbours'),'Neighbours to display',value=30,min=5,max=100),
        shiny::checkboxInput(ns('broad'),'Include broader exploratory GO candidates',FALSE))),
      width=3),
    shiny::mainPanel(shiny::tabsetPanel(id=ns('mode'),
      shiny::tabPanel('Species overview',shiny::uiOutput(ns('overview_summary')),
        shiny::fluidRow(shiny::column(6,plotly::plotlyOutput(ns('overview_qc'),height='340px')),shiny::column(6,plotly::plotlyOutput(ns('overview_coverage'),height='340px'))),
        DT::DTOutput(ns('overview_exclusions'))),
      shiny::tabPanel('What changed?',
        shiny::uiOutput(ns('match_summary')),DT::DTOutput(ns('pairs')),
        shiny::uiOutput(ns('fit_status')),plotly::plotlyOutput(ns('differential'),height='480px'),
        shiny::uiOutput(ns('selected_effect')),plotly::plotlyOutput(ns('replicates'),height='300px'),
        DT::DTOutput(ns('differential_table'))),
      shiny::tabPanel('Which genes move together?',shiny::uiOutput(ns('atlas_note')),
        plotly::plotlyOutput(ns('profiles'),height='350px'),
        plotly::plotlyOutput(ns('network'),height='420px'),DT::DTOutput(ns('neighbours'))),
      shiny::tabPanel('What functions are suggested?',
        shiny::p('Candidate ranks are hypotheses, not probabilities or candidate-level significance. Negative coordination can support a regulatory relationship.'),
        shiny::selectizeInput(ns('pathway'),'Pathway (optional; clears gene filtering)',choices=NULL),
        DT::DTOutput(ns('candidates')),plotly::plotlyOutput(ns('term_performance'),height='400px'),
        DT::DTOutput(ns('members'))),
      selected='Species overview'),width=9)))
}

.translation_cache <- new.env(parent=emptyenv())
translation_cached_read <- function(path) {
  if (!file.exists(path)) stop('Translation data file not found: ',basename(path))
  info<-file.info(path);key<-paste(normalizePath(path),info$size,as.numeric(info$mtime),sep='|')
  if (!exists(key,.translation_cache,inherits=FALSE)) assign(key,readRDS(path),.translation_cache)
  get(key,.translation_cache,inherits=FALSE)
}

translation_server <- function(id, browser_options=NULL, all_experiments=NULL) {
  shiny::moduleServer(id,function(input,output,session) {
    directory<-translation_data_dir(browser_options)
    # Resolve a moving `current` alias once per session, before any lazy reads.
    if(dir.exists(directory))directory<-normalizePath(directory,mustWork=TRUE)
    index<-if(file.exists(file.path(directory,'datasets.csv'))) as.data.frame(data.table::fread(file.path(directory,'datasets.csv'))) else
      if(file.exists(file.path(directory,'atlas.rds')))data.frame(species=translation_cached_read(file.path(directory,if(file.exists(file.path(directory,'dataset.rds')))'dataset.rds' else 'atlas.rds'))$species,directory=directory) else data.frame(species=character(),directory=character())
    if(anyDuplicated(index$species))stop('Translation dataset species identifiers must be unique.')
    if(nrow(index))index$directory<-normalizePath(path.expand(index$directory),mustWork=TRUE)
    shiny::updateSelectInput(session,'species',choices=index$species)
    base<-shiny::reactive({shiny::req(input$species);p<-index$directory[match(input$species,index$species)];shiny::req(length(p)==1L);p})
    bundle<-shiny::reactive({
      rc_access_check()
      b<-translation_cached_read(file.path(base(),if(file.exists(file.path(base(),'atlas_core.rds')))'atlas_core.rds' else 'atlas.rds'))
      if(!identical(b$version,1L))stop('Unsupported Translation dataset version.')
      # Curated atlas summaries must not bypass an access-controlled catalog.
      if(!is.null(rc_access_context())) {
        allowed<-all_experiments$name
        if(!all(b$source_experiments %in% allowed))stop('This atlas includes experiments outside your authorized catalog. Configure an atlas built for that catalog.')
      }
      b
    })
    header<-shiny::reactive({shiny::req(input$species);p<-file.path(base(),'dataset.rds');if(file.exists(p))translation_cached_read(p) else {b<-bundle();list(species=b$species,genes=nrow(b$TE),contexts=ncol(b$TE),source_experiments=b$source_experiments)}})
    check_scope<-function(x){rc_access_check();if(!is.null(rc_access_context())&&!all(x$source_experiments%in%all_experiments$name))stop('This summary includes experiments outside your authorized catalog.')}
    qc_metadata<-shiny::reactive({check_scope(header());p<-file.path(base(),'pair_metadata.rds');if(file.exists(p))translation_cached_read(p) else NULL})
    go_bundle<-shiny::reactive({check_scope(header());p<-file.path(base(),'atlas_go.rds');if(file.exists(p))translation_cached_read(p) else bundle()})
    overview<-shiny::reactive({check_scope(header());p<-file.path(base(),'overview.rds');shiny::validate(shiny::need(file.exists(p),'Prepare overview display sidecars for this release.'));x<-translation_cached_read(p);check_scope(x);x})
    store<-shiny::reactive({bundle();translation_cached_read(file.path(base(),'counts.rds'))})
    mapping<-shiny::reactive({
      p<-if(is.null(input$mapping))file.path(base(),'run_map.csv') else input$mapping$datapath
      as.data.frame(data.table::fread(p,colClasses='character',na.strings=c('','NA')))
    })
    custom_pairs_path<-translation_option(browser_options, 'translation_pairs_file')
    custom_pairs<-!is.null(custom_pairs_path)&&length(custom_pairs_path)&&nzchar(custom_pairs_path)
    pairs<-shiny::reactive({
      bundle()
      if(custom_pairs) {
        if(!is.null(input$mapping))stop('This cohort uses translation_pairs_file. Rebuild that object with your reviewed Run mapping to change pairing.')
        p<-translation_cached_read(path.expand(custom_pairs_path))
        if(!identical(p$species,bundle()$species))stop('Supplied differential cohort and atlas species differ.')
        if(!identical(dimnames(p$Ribo),dimnames(p$RNA)))stop('Supplied paired count identities differ.')
        if(!is.null(rc_access_context()))stop('Custom paired files require an explicitly scoped application; they are disabled in access-controlled sessions.')
        translation_attach_pair_qc(p,qc_metadata())
      } else translation_attach_pair_qc(translation_pair_store(store(),mapping()),qc_metadata())
    })
    gene_initialized<-shiny::reactiveVal(NULL);go_initialized<-shiny::reactiveVal(NULL)
    shiny::observeEvent(list(input$species,input$mode),{
      shiny::req(input$mode,input$species)
      if(input$mode=='Species overview')return()
      if(!identical(gene_initialized(),base())) {
        b<-bundle();labels<-ifelse(is.na(b$genes$symbol)|!nzchar(b$genes$symbol),b$genes$gene_id,paste0(b$genes$symbol,' · ',b$genes$gene_id))
        shiny::updateSelectizeInput(session,'gene',choices=stats::setNames(b$genes$gene_id,labels),selected=if('ATF4'%in%b$genes$symbol)b$genes$gene_id[match('ATF4',b$genes$symbol)] else b$genes$gene_id[1],server=TRUE)
        ctx<-colnames(b$TE);shiny::updateSelectizeInput(session,'contexts',choices=ctx,selected=ctx,server=TRUE)
        gene_initialized(base())
      }
      if(input$mode=='What functions are suggested?'&&!identical(go_initialized(),base())) {
        terms<-go_bundle()$terms;shiny::updateSelectizeInput(session,'pathway',choices=c('Selected gene'='',stats::setNames(terms$GO_ID,paste(terms$GO_name,terms$GO_ID))),server=TRUE);go_initialized(base())
      }
    })
    pair_metadata<-shiny::reactive({
      shiny::req(input$mode=='What changed?')
      if(custom_pairs||!is.null(input$mapping)||is.null(qc_metadata()))return(pairs()$units)
      qc_metadata()$units
    })
    shiny::observeEvent(pair_metadata(),{
      u<-pair_metadata();studies<-sort(unique(u$study[!is.na(u$study)&nzchar(u$study)]))
      has_contrast<-vapply(studies,function(s)length(unique(u$condition[u$study==s & !is.na(u$condition)]))>1L,logical(1))
      shiny::updateSelectInput(session,'study',choices=c('All studies'='__all__',stats::setNames(studies,studies)),selected=if(any(has_contrast))studies[which(has_contrast)[1]] else studies[1])
    })
    study_units<-shiny::reactive({shiny::req(input$study);u<-pair_metadata();if(input$study!='__all__')u<-u[!is.na(u$study)&u$study==input$study,,drop=FALSE];if(!'context'%in%names(u))u$context<-'Unspecified';u})
    shiny::observeEvent(study_units(),{u<-study_units();cs<-sort(unique(u$context));cs<-cs[!is.na(cs)&nzchar(cs)];shiny::updateSelectInput(session,'cohort_context',choices=c('All contexts'='__all__',stats::setNames(cs,cs)),selected=if(length(cs)==1L)cs else '__all__')})
    context_units<-shiny::reactive({u<-study_units();if(!is.null(input$cohort_context)&&input$cohort_context!='__all__')u<-u[!is.na(u$context)&u$context==input$cohort_context,,drop=FALSE];u})
    eligible_units<-shiny::reactive(translation_filter_pair_qc(context_units(),minimum=if(is.null(input$pair_r2_min)).188 else input$pair_r2_min,enabled=isTRUE(input$pair_qc),include_unscored=isTRUE(input$include_unscored)))
    shiny::observeEvent(list(eligible_units(),input$all_units),{
      u<-eligible_units();if(isTRUE(input$all_units))return()
      selected<-shiny::isolate(input$units);selected<-intersect(selected,u$biosample_id)
      shiny::updateSelectizeInput(session,'units',choices=stats::setNames(u$biosample_id,paste(u$biosample_id,u$condition,u$context,sep=' · ')),selected=selected,server=TRUE)
    })
    selected_units<-shiny::reactive({u<-eligible_units();if(isTRUE(input$all_units))return(u);shiny::req(input$units);u[u$biosample_id %in% input$units,,drop=FALSE]})
    output$unit_count<-shiny::renderUI({u<-eligible_units();shiny::p(paste(nrow(u),'eligible of',nrow(context_units()),'profiles after pairing QC. Biological replicates stay separate.'))})
    shiny::observeEvent(eligible_units(),{
      cs<-sort(unique(eligible_units()$condition));cs<-cs[!is.na(cs)&nzchar(cs)]
      if(!length(cs)){shiny::updateSelectInput(session,'reference',choices=character(),selected=character());shiny::updateSelectInput(session,'comparison',choices=character(),selected=character());return()}
      ref<-if('WT'%in%cs)'WT' else cs[1]
      shiny::updateSelectInput(session,'reference',choices=cs,selected=ref)
      shiny::updateSelectInput(session,'comparison',choices=cs,selected=head(setdiff(cs,ref),1))
    })
    output$status<-shiny::renderUI({
      if(!nrow(index))return(shiny::p('No Translation dataset configured. Set translation_data_dir or RIBOCRYPT_TRANSLATION_DIR.'))
      b<-header();check_scope(b);shiny::p(paste(b$genes,'genes ·',b$contexts,'contexts'),shiny::br(),paste('Species:',b$species))
    })
    output$overview_summary<-shiny::renderUI({
      shiny::req(input$mode=='Species overview');x<-overview();labels<-c(curated_pairs='Curated assay pairs',count_ready_pairs='Pairs with counts',atlas_pairs='Atlas pairs',contexts='Contexts',studies='Atlas studies')
      cards<-lapply(names(labels),function(key)shiny::div(class='translation-overview-card',shiny::strong(format(x$counts[[key]],big.mark=',')),labels[[key]]))
      shiny::tagList(shiny::h4(paste(x$species,'matched RNA/Ribo data')),shiny::div(class='translation-overview-cards',cards),
        shiny::p('Curated correspondence is metadata-supported; exact biological material is not established for every pair. Counts-ready pairs pass additional pairing, CDS and context-support gates before atlas inclusion.'),
        shiny::p(paste('Reference atlas:',x$counts[['atlas_pairs']],'pairs; minimum CLR-count R²',x$reference_cutoff,'. Differential QC filters narrow this available cohort; they do not rebuild the saved GO reference.')))
    })
    output$overview_qc<-plotly::renderPlotly({
      shiny::req(input$mode=='Species overview');x<-overview();z<-x$pairs;z<-z[is.finite(z$matched_R2),];z$cohort<-ifelse(z$included_in_primary,'In reference atlas','Outside reference atlas')
      plotly::plot_ly(z,x=~matched_R2,color=~cohort,type='histogram',nbinsx=35) |>
        plotly::layout(barmode='overlay',xaxis=list(title='Matched RNA–Ribo R² (CLR counts)'),yaxis=list(title='Pairs'),shapes=list(list(type='line',x0=x$reference_cutoff,x1=x$reference_cutoff,y0=0,y1=1,yref='paper',line=list(dash='dash',color='gray'))))
    })
    output$overview_coverage<-plotly::renderPlotly({
      shiny::req(input$mode=='Species overview');x<-overview();z<-head(x$contexts[order(-x$contexts$n_pairs),],12)
      plotly::plot_ly(z,x=~n_pairs,y=~profile_context,type='bar',orientation='h',text=~paste(n_studies,'studies')) |>
        plotly::layout(xaxis=list(title='Atlas pairs'),yaxis=list(title='Top contexts',autorange='reversed'),margin=list(l=140))
    })
    output$overview_exclusions<-DT::renderDT(DT::datatable({shiny::req(input$mode=='Species overview');overview()$exclusions},rownames=FALSE,options=list(pageLength=6,dom='tip')),server=TRUE)
    output$match_summary<-shiny::renderUI({
      u<-selected_units();shiny::p(paste(nrow(u),'matched profiles selected for this contrast. Counts load only when a fit is requested.'),shiny::br(),
        if(custom_pairs)'This cohort comes from your paired-count object. Confirm independent biological replication and material identity before inference.' else 'The supplied atlas mapping is metadata-supported; shared biological material has not been established. Confirm replication and material identity before inferential paired analysis.')
    })
    output$pairs<-DT::renderDT({u<-selected_units();DT::datatable(u,rownames=FALSE,options=list(pageLength=5,scrollX=TRUE))},server=TRUE)
    active_key<-shiny::reactiveVal(NULL)
    result<-shiny::reactiveVal(NULL);fit_message<-shiny::reactiveVal('Select a study and run a descriptive or replicate-aware differential analysis.')
    selected_gene<-shiny::reactiveVal(NULL)
    task<-NULL
    if(nrow(index)>0L && requireNamespace('future',quietly=TRUE) && requireNamespace('promises',quietly=TRUE) && utils::packageVersion('shiny')>='1.8.1') {
      fit_function<-translation_fit
      task<-shiny::ExtendedTask$new(function(payload) {
        promises::future_promise({fit_function(payload)},globals=list(fit_function=fit_function,payload=payload),packages='DESeq2',seed=TRUE)
      })
      shiny::observeEvent(task$status(),{
        status<-task$status()
        if(status=='running')fit_message('Fitting in a background worker; atlas browsing remains available.')
        if(status=='success'){
          ans<-task$result();if(!identical(ans$source_key,active_key()))return();result(ans);fit_message(paste('Finished:',ans$retained_features,'features.',ans$model))
          # Disk result cache is isolated by the app's configured cohort release.
          if(isTRUE(ans$cacheable))try(saveRDS(ans,file.path(base(),'fit_cache',paste0(ans$source_key,'.rds'))),silent=TRUE)
        }
        if(status=='error')tryCatch(task$result(),error=function(e)fit_message(conditionMessage(e)))
      })
    }
    shiny::observeEvent(input$fit,{
      rc_access_check();result(NULL)
      tryCatch({
        u<-selected_units();u<-u[u$condition %in% c(input$reference,input$comparison),,drop=FALSE]
        if(!nrow(u))stop('No biological units in this contrast.')
        if(length(input$reference)!=1L||length(input$comparison)!=1L||identical(input$reference,input$comparison))stop('Select distinct reference and comparison conditions.')
        if(!isTRUE(input$descriptive)&&(!isTRUE(input$replication)||any(table(factor(u$condition,levels=c(input$reference,input$comparison)))<2L)))stop('Confirm at least two independent biological units per condition before inference.')
        if(!isTRUE(input$descriptive) && length(unique(u$context))>1L)stop('Select a single biological context before differential inference.')
        p<-pairs()
        payload<-list(Ribo=p$Ribo[,u$biosample_id,drop=FALSE],RNA=p$RNA[,u$biosample_id,drop=FALSE],units=u,
                      reference=input$reference,comparison=input$comparison,descriptive=input$descriptive,cacheable=!custom_pairs&&is.null(input$mapping)&&is.null(rc_access_context()),
                      reviewed_replication=input$replication,paired=input$paired,batch=if(length(unique(u$study))>1L)'study' else NULL,genes=p$genes,count_region=p$count_region)
        if(payload$paired && !isTRUE(input$replication))stop('Review independent biological replication before requesting an exact-material paired model.')
        if(payload$paired && !custom_pairs && is.null(input$mapping))stop('Upload a reviewed same-material mapping; the default atlas does not establish specimen identity.')
        payload$source_key<-digest::digest(list(payload,source=bundle()$source_run,reference=p$reference,engine=body(translation_fit)),algo='sha256')
        active_key(payload$source_key);selected_gene(NULL)
        cache<-file.path(base(),'fit_cache');dir.create(cache,showWarnings=FALSE)
        path<-file.path(cache,paste0(payload$source_key,'.rds'))
        if(isTRUE(payload$cacheable) && file.exists(path)){ans<-readRDS(path);result(ans);fit_message('Loaded cached differential analysis.');return()}
        if(is.null(task))stop('Install future and promises, and use Shiny >= 1.8.1 for background fits.')
        if(inherits(future::plan(),'sequential'))future::plan(future::multisession,workers=I(1))
        task$invoke(payload)
      },error=function(e)fit_message(conditionMessage(e)))
    })
    shiny::observeEvent(list(input$species,input$mapping,input$pair_qc,input$pair_r2_min,input$include_unscored),{active_key(NULL);result(NULL)},ignoreInit=TRUE)
    output$fit_status<-shiny::renderUI(shiny::p(fit_message()))
    filtered<-shiny::reactive({
      rc_access_check();r<-result();shiny::req(r);z<-r$table
      z$Evidence<-if(r$descriptive)ifelse(abs(z$deltaTE)>=input$effect,'Effect filter passed','Below effect filter') else
        ifelse(!is.na(z$FDR)&z$FDR<=input$fdr&abs(z$deltaTE)>=input$effect,'Differential TE','Other / insufficient evidence')
      z
    })
    output$differential<-plotly::renderPlotly({
      z<-filtered();r<-result();p<-plotly::plot_ly(z,x=~RNA_log2FC,y=~Ribo_log2FC,key=~gene_id,type='scattergl',mode='markers',color=~Evidence,
        text=~paste(symbol,gene_id,'ΔTE:',signif(deltaTE,3),'FDR:',signif(FDR,3)),hoverinfo='text',source=session$ns('dte'))
      plotly::event_register(plotly::layout(p,xaxis=list(title=paste0('RNA log₂ fold change (',r$comparison,' / ',r$reference,')')),
        yaxis=list(title=paste0('Ribo log₂ fold change (',r$comparison,' / ',r$reference,')')),
        shapes=list(list(type='line',xref='x',yref='y',x0=-15,y0=-15,x1=15,y1=15,line=list(dash='dash',color='gray')))),'plotly_click')
    })
    click_event<-shiny::reactive({shiny::req(result());plotly::event_data('plotly_click',source=session$ns('dte'))})
    shiny::observeEvent(click_event(),{e<-click_event();if(!is.null(e$key))selected_gene(e$key[1])})
    output$differential_table<-DT::renderDT(DT::datatable(filtered(),rownames=FALSE,selection='single',options=list(pageLength=10,scrollX=TRUE)),server=TRUE)
    shiny::observeEvent(input$differential_table_rows_selected,{i<-input$differential_table_rows_selected;if(length(i))selected_gene(filtered()$gene_id[i[1]])})
    current_gene<-shiny::reactive({if(!is.null(selected_gene()) && input$mode=='What changed?')selected_gene() else input$gene})
    output$selected_effect<-shiny::renderUI({
      r<-result();g<-current_gene();shiny::req(r,g);z<-r$table[r$table$gene_id==g,,drop=FALSE];shiny::req(nrow(z)==1)
      # Use a validated study-specific experiment when it is available.
      u<-r$units;exp<-if('ribo_experiment'%in%names(u))u$ribo_experiment[1] else ''
      link<-if(nzchar(exp))make_rc_url(gene_id=g,exp=exp,host=getHostFromURL(session),plot_on_start=TRUE) else NULL
      shiny::tagList(shiny::strong(paste(z$symbol,g)),shiny::p(paste('ΔTE:',signif(z$deltaTE,3),'95% interval:',signif(z$CI_low,3),'to',signif(z$CI_high,3),'FDR:',signif(z$FDR,3))),
        shiny::p(r$normalization,paste('Regions:',r$regions)),if(!is.null(link))shiny::a('Open footprints in RiboCrypt',href=link,target='_blank'))
    })
    output$replicates<-plotly::renderPlotly({
      r<-result();g<-current_gene();shiny::req(r,g,g%in%rownames(r$normalized_RNA));u<-r$units
      z<-data.frame(biosample_id=rep(u$biosample_id,2),condition=rep(u$condition,2),assay=rep(c('RNA','Ribo'),each=nrow(u)),
        value=log2(c(r$normalized_RNA[g,],r$normalized_Ribo[g,])+.5))
      plotly::plot_ly(z,x=~condition,y=~value,color=~assay,text=~biosample_id,type='scatter',mode='markers') |>
        plotly::layout(yaxis=list(title='Assay-normalized abundance, log₂(count + 0.5)'))
    })
    output$download<-shiny::downloadHandler(filename=function()'differential_TE.csv',content=function(file){rc_access_check(capability='download');data.table::fwrite(filtered(),file)})
    neighbours<-shiny::reactive({shiny::req(input$mode=='Which genes move together?');rc_access_check();shiny::req(input$gene,input$contexts);translation_neighbours(bundle(),input$gene,input$contexts,max(5L,min(100L,as.integer(input$n_neighbours))))})
    output$atlas_note<-shiny::renderUI({shiny::req(input$mode=='Which genes move together?');b<-bundle();shiny::tagList(shiny::p('Explore coordinated and opposing translation patterns across cell types and tissues. Neighbours reflect the selected contexts; GO predictions use the full reference atlas.'),shiny::tags$details(shiny::tags$summary('Method and data source'),shiny::p(b$TE_definition,shiny::br(),b$pooling,shiny::br(),paste('Source release:',basename(b$source_run)))))})
    output$neighbours<-DT::renderDT(DT::formatRound(DT::datatable(neighbours(),rownames=FALSE,options=list(pageLength=10)),c('TE_rho','RNA_rho'),digits=3),server=TRUE)
    output$profiles<-plotly::renderPlotly({
      shiny::req(input$mode=='Which genes move together?');b<-bundle();n<-neighbours();ids<-c(input$gene,head(n$gene_id,5));x<-b$TE[ids,input$contexts,drop=FALSE]
      x<-x-rowMeans(x)
      z<-data.frame(gene_id=rep(rownames(x),ncol(x)),context=rep(colnames(x),each=nrow(x)),value=as.vector(x))
      z$gene<-b$genes$symbol[match(z$gene_id,b$genes$gene_id)];missing<-is.na(z$gene)|!nzchar(z$gene);z$gene[missing]<-z$gene_id[missing]
      plotly::plot_ly(z,x=~context,y=~value,color=~gene,type='scatter',mode='lines+markers') |>
        plotly::layout(yaxis=list(title='Compositional TE, centered per gene'),xaxis=list(title='Context',tickangle=-45))
    })
    output$network<-plotly::renderPlotly({
      n<-neighbours();shiny::req(nrow(n)>0);angles<-seq(0,2*pi,length.out=nrow(n)+1)[-1]
      x<-cos(angles);y<-sin(angles)
      edges<-data.frame(x=as.vector(rbind(0,x,NA)),y=as.vector(rbind(0,y,NA)))
      p<-plotly::plot_ly(edges,x=~x,y=~y,type='scatter',mode='lines',hoverinfo='skip',line=list(color='lightgray'),showlegend=FALSE)
      nodes<-data.frame(x=c(0,x),y=c(0,y),label=c(input$gene,ifelse(is.na(n$symbol),n$gene_id,n$symbol)),rho=c(1,n$TE_rho),gene_id=c(input$gene,n$gene_id))
      plotly::add_trace(p,data=nodes,x=~x,y=~y,type='scatter',mode='markers',inherit=FALSE,color=~rho,colors=c('#3b6fb6','#f6f6f6','#d84a3a'),text=~paste(label,'signed TE rho:',signif(rho,3)),hoverinfo='text',showlegend=FALSE) |>
        plotly::layout(xaxis=list(visible=FALSE),yaxis=list(visible=FALSE,scaleanchor='x'))
    })
    candidates<-shiny::reactive({
      shiny::req(input$mode=='What functions are suggested?');rc_access_check();b<-go_bundle();z<-if(isTRUE(input$broad))b$candidates else b$retained
      if(!is.null(input$pathway)&&nzchar(input$pathway))z<-z[z$GO_ID==input$pathway,] else {shiny::req(input$gene);z<-z[z$gene_id==input$gene,]}
      z
    })
    output$candidates<-DT::renderDT({z<-as.data.frame(candidates());cols<-intersect(c('symbol','gene_id','GO_name','GO_ID','TE_rank','TE_percentile','RNA_percentile','term_TE_AUROC','term_RNA_AUROC','TE_mean_abs_rho','all_models_pass'),names(z));DT::datatable(z[,cols,drop=FALSE],rownames=FALSE,selection='single',options=list(pageLength=10,scrollX=TRUE))},server=TRUE)
    shiny::observeEvent(input$candidates_rows_selected,{i<-input$candidates_rows_selected;if(length(i)){g<-candidates()$gene_id[i[1]];shiny::updateSelectizeInput(session,'gene',selected=g)}})
    output$term_performance<-plotly::renderPlotly({shiny::req(input$mode=='What functions are suggested?');b<-go_bundle();z<-b$terms;plotly::plot_ly(z,x=~RNA_AUROC,y=~TE_AUROC,text=~paste(GO_name,GO_ID,'members:',n_genes),type='scattergl',mode='markers') |>
      plotly::layout(xaxis=list(title='RNA leave-one-gene-out AUROC'),yaxis=list(title='TE leave-one-gene-out AUROC'),shapes=list(list(type='line',x0=0,y0=0,x1=1,y1=1,line=list(dash='dash',color='gray'))))})
    output$members<-DT::renderDT({
      shiny::req(input$mode=='What functions are suggested?');b<-go_bundle();term<-input$pathway
      if(is.null(term)||!nzchar(term)){z<-candidates();shiny::validate(shiny::need(nrow(z)>0,'No candidate rows for this gene in the selected evidence set. Choose a pathway or include exploratory candidates.'));term<-z$GO_ID[1]}
      m<-b$memberships[b$memberships$GO_ID==term,];m$status<-'Known GO member';DT::datatable(m,rownames=FALSE,options=list(pageLength=10,scrollX=TRUE))
    },server=TRUE)
  })
}

#' Preview the Translation analysis page
#' @param data_dir Directory produced by prepare_translation_atlas.R, or a
#'   directory containing datasets.csv (species,directory) for multiple releases.
#' @param options Shiny application options.
#' @param paired_counts_file Optional RDS from ribocrypt_translation_pairs for
#'   differential analysis of custom Run-resolved collections of the same species.
#' @return A Shiny application; no public server is modified.
#' @export
ribocrypt_translation_app <- function(data_dir=Sys.getenv('RIBOCRYPT_TRANSLATION_DIR'),options=list(), paired_counts_file=NULL) {
  shiny::shinyApp(shiny::fluidPage(shiny::titlePanel('RiboCrypt · Translation'),
    shiny::tabsetPanel(translation_ui('translation'))),
    function(input,output,session)translation_server('translation',c(translation_data_dir=data_dir,translation_pairs_file=if(is.null(paired_counts_file))'' else paired_counts_file)),options=options)
}
