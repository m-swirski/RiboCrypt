translation_fixture <- function() {
  ids<-paste0('g',1:60);libs<-c('r1','r2','a1','a2')
  counts<-matrix(20L,60,2,dimnames=list(ids,libs[1:2]))
  A<-counts;colnames(A)<-libs[3:4]
  store<-list(Ribo=counts,RNA=A,
    run_members=data.frame(Run=c('SRR1','SRR2','SRR3','SRR4','SRR5'),library_id=c('r1','r1','a1','r2','a2'),assay=c('RIBO','RIBO','RNA','RIBO','RNA')),
    units=data.frame(biosample_id=c('b1','b2'),ribo_library_id=c('r1','r2'),rna_library_id=c('a1','a2'),condition=c('WT','treated'),study='s',replicate=c('1','2'),match_level='reviewed'),
    genes=data.frame(gene_id=ids,symbol=ids),count_region='CDS/mRNA',reference='ref')
  map<-data.frame(Run=c('SRR1','SRR2','SRR3','SRR4','SRR5'),biosample_id=c('b1','b1','b1','b2','b2'))
  list(store=store,map=map)
}

test_that('technical components and keyed ordering preserve biological units',{
  f<-translation_fixture();p<-translation_pair_store(f$store,f$map)
  expect_identical(colnames(p$Ribo),c('b1','b2'));expect_equal(nrow(p$units),2)
  expect_equal(p$inventory$technical_runs[p$inventory$library_id=='r1'],2)
  f$store$Ribo<-f$store$Ribo[,2:1];f$store$RNA<-f$store$RNA[,2:1]
  q<-translation_pair_store(f$store,f$map[5:1,]);expect_identical(q$Ribo,p$Ribo);expect_identical(q$RNA,p$RNA)
  f$map$biosample_id[1]<-'conflict';expect_error(translation_pair_store(f$store,f$map),'conflicting')
})

test_that('NA stays unmatched and incomplete or duplicate mappings cannot inflate n',{
  f<-translation_fixture();f$map$biosample_id[5]<-NA
  p<-translation_pair_store(f$store,f$map);expect_equal(nrow(p$units),1);expect_equal(p$unmatched_biosamples,'b2')
  f<-translation_fixture();expect_error(translation_pair_store(f$store,f$map[-2,]),'Incomplete')
  expect_error(translation_pair_store(f$store,rbind(f$map,f$map[1,])),'unique')
  f$map$biosample_id<-'one';expect_error(translation_pair_store(f$store,f$map),'Multiple independent')
  f<-translation_fixture();f$map$Run[1]<-'unknown';expect_error(translation_pair_store(f$store,f$map),'absent')
})

test_that('signed atlas query agrees with literal proportionality and rejects missing contexts',{
  x<-rbind(g1=c(1,2,4,5),g2=c(2,4,8,10),g3=c(-1,-2,-4,-5),g4=rep(0,4))
  colnames(x)<-paste0('c',1:4);b<-list(TE=x,RNA=x,genes=data.frame(gene_id=rownames(x),symbol=rownames(x)))
  z<-translation_neighbours(b,'g1',n=10)
  expect_equal(z$gene_id,c('g3','g2'));expect_equal(z$TE_rho,c(-1,.8))
  expect_false('g4'%in%z$gene_id);expect_error(translation_neighbours(b,'g1',contexts='c1'),'three')
  expect_error(translation_neighbours(b,'unknown'),'not in')
})

translation_fit_fixture <- function() {
  set.seed(7);G<-250;N<-6
  baseline<-rep(200,G);A<-matrix(rnbinom(G*N,mu=baseline,size=30),G,N)
  R<-matrix(rnbinom(G*N,mu=baseline,size=30),G,N);R[1:20,4:6]<-R[1:20,4:6]*4L
  dimnames(A)<-dimnames(R)<-list(paste0('g',1:G),paste0('b',1:N))
  list(Ribo=R,RNA=A,units=data.frame(biosample_id=colnames(A),condition=rep(c('WT','treated'),each=3)),
       genes=data.frame(gene_id=rownames(A),symbol=rownames(A)),reference='WT',comparison='treated',
       descriptive=FALSE,reviewed_replication=TRUE,paired=FALSE,count_region='CDS/mRNA',source_key='test')
}

test_that('inference guards reject fractional counts, missing replication and confounding',{
  skip_if_not_installed('DESeq2');p<-translation_fit_fixture()
  p$Ribo[1,1]<-.5;expect_error(translation_fit(p),'integer');p<-translation_fit_fixture()
  p$reviewed_replication<-FALSE;expect_error(translation_fit(p),'replication')
  p<-translation_fit_fixture();p$batch<-'batch';p$units$batch<-p$units$condition
  expect_error(translation_fit(p),'confounded')
  p<-translation_fit_fixture();p$comparison<-NULL;expect_error(translation_fit(p),'distinct')
})

test_that('differential interaction detects a Ribo-specific spike and reverses with contrast',{
  skip_if_not_installed('DESeq2');p<-translation_fit_fixture();r<-translation_fit(p)
  expect_gt(median(r$table$deltaTE[1:20]),1.5);expect_lt(median(r$table$FDR[1:20]),.05)
  q<-p;q$reference<-'treated';q$comparison<-'WT';s<-translation_fit(q)
  expect_equal(s$table$deltaTE,-r$table$deltaTE,tolerance=1e-5)
  p$paired<-TRUE;z<-translation_fit(p);expect_gt(median(z$table$deltaTE[1:20]),1.5)
  p$descriptive<-TRUE;d<-translation_fit(p);expect_true(all(is.na(d$table$FDR)))
})

test_that('Translation UI is explicit about contrast and is wired to lazy startup',{
  html<-as.character(translation_ui('te'))
  expect_match(html,'Reference (Group 1)',fixed=TRUE);expect_match(html,'Which genes move together?',fixed=TRUE)
  expect_match(paste(deparse(body(analysis_ui)),collapse=' '),'translation_ui',fixed=TRUE)
  expect_match(paste(deparse(body(RiboCrypt_app)),collapse=' '),'Translation',fixed=TRUE)
})

test_that('collection adapter sums reviewed sequencing lanes and joins regional features by ID',{
  r<-list(Run=paste0('R',1:4),condition=c('WT','WT','treated','treated'),rep=c('1','1','2','2'))
  a<-list(Run=paste0('A',1:4),condition=r$condition,rep=r$rep)
  R<-matrix(20L,60,4,dimnames=list(paste0('g',1:60),r$Run));A<-rbind(R,ncRNA=rep(10L,4));colnames(A)<-a$Run
  map<-data.frame(Run=c(r$Run,a$Run),biosample_id=rep(c('b1','b1','b2','b2'),2),assay_library_id=c('r1','r1','r2','r2','a1','a1','a2','a2'))
  testthat::local_mocked_bindings(organism=function(x)'Homo sapiens',.package='BiocGenerics')
  testthat::with_mocked_bindings({
    p<-ribocrypt_translation_pairs(r,a,map,R,A[,4:1])
    expect_equal(nrow(p$units),2L);expect_identical(dimnames(p$Ribo),dimnames(p$RNA))
    expect_true(all(p$Ribo==40));expect_true(all(p$RNA==40));expect_identical(p$species,'Homo sapiens')
  },findFa=function(x)'reference',name=function(x)'study',.package='ORFik')
})

test_that('precomputed neighbours match live queries without data.table name shadowing',{
  x<-rbind(g1=c(1,2,4,5),g2=c(2,4,8,10),g3=c(-1,-2,-4,-5));colnames(x)<-paste0('c',1:4)
  b<-list(TE=x,RNA=x,genes=data.frame(gene_id=rownames(x),symbol=rownames(x)))
  live<-translation_neighbours(b,'g1',n=2)
  pre<-data.table::data.table(query='g1',gene_id=live$gene_id,TE_rho=live$TE_rho)
  expect_identical(translation_neighbours(b,'g1',n=2,precomputed=pre),live)
})

test_that('optional Translation settings accept the main launcher character vector',{
  withr::local_envvar(c(RIBOCRYPT_TRANSLATION_DIR='fallback-atlas'))
  old<-c(default_experiment='all_merged-Homo_sapiens',plot_on_start='FALSE')
  expect_identical(translation_data_dir(old),'fallback-atlas')
  expect_identical(translation_data_dir(c()),'fallback-atlas')
  expect_identical(translation_data_dir(NULL),'fallback-atlas')
  expect_identical(translation_data_dir(c(translation_data_dir=NA_character_)),'fallback-atlas')
  expect_identical(translation_data_dir(c(old,translation_data_dir='configured-atlas')),'configured-atlas')
  expect_identical(translation_option(old,'translation_pairs_file'),'')
  expect_identical(translation_option(list(translation_pairs_file=NULL),'translation_pairs_file'),'')
})

test_that('an unconfigured Translation module does not break another analysis tab',{
  withr::local_envvar(c(RIBOCRYPT_TRANSLATION_DIR=''))
  shiny::testServer(translation_server,args=list(browser_options=c(default_experiment='example')),
    {expect_match(output$status$html,'No Translation dataset configured',fixed=TRUE)})
})

test_that('pair QC is calculated on keyed CLR components and missing scores stay missing',{
  R<-matrix(c(20,50,100,80,5,200,2,60),4,2,dimnames=list(paste0('g',1:4),c('b1','b2')))
  A<-matrix(c(30,70,90,100,12,150,6,70),4,2,dimnames=dimnames(R));p<-list(Ribo=R,RNA=A,reference='ref')
  z<-translation_pair_r2(p)
  expect_equal(z$R2,vapply(1:2,function(j)cor(log(pmax(R[,j],1)),log(pmax(A[,j],1)))^2,numeric(1)))
  qc<-list(reference='ref',feature_ids=c('g1','g2'),other_ids=c('g3','g4'),definition='frozen components')
  z<-translation_pair_r2(p,qc)
  expect_equal(z$R2,vapply(1:2,function(j)cor(log(c(R[1:2,j],sum(R[3:4,j]))),log(c(A[1:2,j],sum(A[3:4,j]))))^2,numeric(1)))
  u<-data.frame(biosample_id=c('a','b','c'),matched_R2=c(.187,.188,NA))
  expect_identical(translation_filter_pair_qc(u)$biosample_id,'b')
  expect_identical(translation_filter_pair_qc(u,include_unscored=TRUE)$biosample_id,c('b','c'))
  expect_identical(translation_filter_pair_qc(u,enabled=FALSE),u)
  expect_error(translation_filter_pair_qc(u,minimum=2),'between')
})

test_that('remapping an RNA library cannot inherit the old pairing QC',{
  f<-translation_fixture();p<-translation_pair_store(f$store,f$map)
  p$RNA[,1]<-c(rep(10,30),rep(300,30));p$Ribo[,1]<-c(rep(15,30),rep(400,30))
  p$RNA[,2]<-c(rep(300,30),rep(10,30));p$Ribo[,2]<-c(rep(400,30),rep(15,30))
  q<-list(reference='other-reference',units=transform(f$store$units,matched_R2=c(.7,.8),matching_score=c(.2,.3),Ribo_CDS_fraction_proxy=1,qc_rna_library_id=c('a1','wrong')))
  z<-translation_attach_pair_qc(p,q)
  expect_equal(z$units$matched_R2[1],.7)
  expect_gt(z$units$matched_R2[2],.99)
})

test_that('overview startup reads no atlas matrices, counts, neighbour graph or worker',{
  dir<-tempfile();dir.create(dir);withr::defer(unlink(dir,recursive=TRUE))
  saveRDS(list(version=1L,species='Homo sapiens',genes=1,contexts=1,source_experiments=character()),file.path(dir,'dataset.rds'))
  saveRDS(list(),file.path(dir,'atlas.rds')) # legacy marker: must not be deserialized
  saveRDS(list(species='Homo sapiens',counts=c(curated_pairs=2,count_ready_pairs=1,atlas_pairs=1,contexts=1,studies=1),source_experiments=character(),reference_cutoff=.188),file.path(dir,'overview.rds'))
  rm(list=ls(.translation_cache),envir=.translation_cache)
  shiny::testServer(translation_server,args=list(browser_options=c(translation_data_dir=dir)),{
    session$setInputs(species='Homo sapiens',mode='Species overview');session$flushReact()
    expect_match(output$status$html,'1 genes',fixed=TRUE)
    expect_match(output$overview_summary$html,'Curated assay pairs',fixed=TRUE)
    expect_false(any(grepl('counts.rds|atlas.rds|atlas_core.rds|top_neighbours.rds',ls(.translation_cache))))
  })
})
