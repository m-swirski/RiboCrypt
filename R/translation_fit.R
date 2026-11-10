# Self-contained so this function can execute in a background R worker while
# the parent application is loaded with devtools rather than installed.
translation_fit <- function(payload) {
  R <- payload$Ribo; A <- payload$RNA; units <- payload$units
  if (!identical(dimnames(R),dimnames(A)) || !identical(colnames(R),as.character(units$biosample_id))) stop('Counts and biological units are not aligned.')
  if (anyDuplicated(units$biosample_id) || any(units$condition %in% NA)) stop('Biological identities/conditions are missing or duplicated.')
  if (!all(is.finite(R)) || !all(is.finite(A)) || any(R<0) || any(A<0)) stop('Counts must be finite and nonnegative.')
  if (any(R!=round(R)) || any(A!=round(A))) stop('Count inference requires integer counts. Fractional assignments must use a validated alternative; do not silently round.')
  ref<-payload$reference;cmp<-payload$comparison
  if (length(ref)!=1L || length(cmp)!=1L || is.na(ref) || is.na(cmp) || !nzchar(ref) || !nzchar(cmp) || identical(ref,cmp) || !all(units$condition %in% c(ref,cmp))) stop('Select distinct reference and comparison conditions.')
  condition<-factor(units$condition,levels=c(ref,cmp)); n<-table(condition)
  descriptive<-isTRUE(payload$descriptive)
  if (!descriptive && (!isTRUE(payload$reviewed_replication) || any(n<2L))) stop('Confirm independent biological replication; inference needs at least two units per condition.')
  if (!requireNamespace('DESeq2',quietly=TRUE)) stop('Install DESeq2 for differential Translation analysis.')
  enough<-max(1L,min(n));keep<-rowSums(A>=10)>=enough & rowSums(R>=10)>=enough
  if (sum(keep)<20L) stop('Fewer than 20 features have supported counts in both assays.')
  R<-R[keep,,drop=FALSE]; A<-A[keep,,drop=FALSE]
  # Feature filtering belongs to the fit key, never the gene search or display.
  factors<-function(x) {
    gm<-exp(rowSums(log(ifelse(x>0,x,1)))/ncol(x));gm[rowSums(x)==0]<-0
    sf<-DESeq2::estimateSizeFactorsForMatrix(x,geoMeans=gm)
    if (any(!is.finite(sf)) || any(sf<=0)) stop('Cannot estimate assay depth factors.')
    sf/exp(mean(log(sf)))
  }
  sfA<-factors(A);sfR<-factors(R)
  normA<-sweep(A,2,sfA,'/');normR<-sweep(R,2,sfR,'/')
  avg<-function(x,g)rowMeans(x[,condition==g,drop=FALSE])
  lfcA<-log2((avg(normA,cmp)+.5)/(avg(normA,ref)+.5));lfcR<-log2((avg(normR,cmp)+.5)/(avg(normR,ref)+.5))
  ans<-data.frame(gene_id=rownames(A),RNA_log2FC=lfcA,Ribo_log2FC=lfcR,
                  deltaTE=lfcR-lfcA,SE=NA_real_,pvalue=NA_real_,FDR=NA_real_,
                  RNA_FDR=NA_real_,Ribo_FDR=NA_real_,mean_RNA=rowMeans(normA),mean_Ribo=rowMeans(normR))
  if (!descriptive) {
    cd<-data.frame(condition=condition,row.names=colnames(A)); batch<-payload$batch
    main_design<-~condition
    if (!is.null(batch) && nzchar(batch)) {
      if (!batch %in% names(units) || anyNA(units[[batch]])) stop('Selected batch metadata is missing.')
      cd$batch<-factor(units[[batch]])
      if(nlevels(cd$batch)>1L)main_design<-~batch+condition
    }
    check<-function(form,data) {
      mm<-stats::model.matrix(form,data)
      if (qr(mm)$rank<ncol(mm) || nrow(mm)<=ncol(mm)) stop('Design is confounded or has no residual degrees of freedom.')
    }
    check(main_design,cd)
    fit<-function(counts,data,design,sf) {
      z<-DESeq2::DESeqDataSetFromMatrix(counts,colData=data,design=design)
      DESeq2::sizeFactors(z)<-sf
      DESeq2::DESeq(z,quiet=TRUE,fitType='local',minReplicatesForReplace=Inf)
    }
    dA<-fit(A,cd,main_design,sfA);dR<-fit(R,cd,main_design,sfR)
    a<-DESeq2::results(dA,contrast=c('condition',cmp,ref));r<-DESeq2::results(dR,contrast=c('condition',cmp,ref))
    joint<-cbind(A,R);colnames(joint)<-paste0(rep(c('RNA::','RIBO::'),each=ncol(A)),rep(colnames(A),2))
    jc<-cd[rep(seq_len(nrow(cd)),2),,drop=FALSE];rownames(jc)<-colnames(joint)
    jc$assay<-factor(rep(c('RNA','RIBO'),each=ncol(A)),levels=c('RNA','RIBO'))
    jc$dte<-as.numeric(jc$assay=='RIBO' & jc$condition==cmp)
    jd<-if('batch'%in%names(jc) && nlevels(jc$batch)>1L) ~batch+condition+assay+dte else ~condition+assay+dte
    if (isTRUE(payload$paired)) {
      jc$sample<-factor(rep(units$biosample_id,2));jd<-~sample+assay+dte
    }
    check(jd,jc)
    j<-fit(joint,jc,jd,c(sfA,sfR));te<-DESeq2::results(j,name='dte')
    # All displayed model effects are unshrunk and have their own provenance.
    ans$RNA_log2FC<-a$log2FoldChange;ans$Ribo_log2FC<-r$log2FoldChange
    ans$deltaTE<-te$log2FoldChange;ans$SE<-te$lfcSE;ans$pvalue<-te$pvalue;ans$FDR<-te$padj
    ans$RNA_FDR<-a$padj;ans$Ribo_FDR<-r$padj
  }
  genes<-payload$genes;ans$symbol<-genes$symbol[match(ans$gene_id,genes$gene_id)]
  ans$CI_low<-ans$deltaTE-1.96*ans$SE;ans$CI_high<-ans$deltaTE+1.96*ans$SE
  list(table=ans,normalized_RNA=normA,normalized_Ribo=normR,units=units,
       reference=ref,comparison=cmp,descriptive=descriptive,paired=isTRUE(payload$paired),
       normalization='Separate assay median-ratio size factors (poscounts), geometric mean one',
       regions=payload$count_region,retained_features=sum(keep),excluded_features=sum(!keep),
       model=if(descriptive)'Descriptive log ratios with 0.5 normalized-count pseudocount; no significance' else if(isTRUE(payload$paired))'sample + assay + treatment:Ribo indicator' else 'condition + assay + treatment:Ribo indicator; optional batch',
       source_key=payload$source_key,cacheable=isTRUE(payload$cacheable))
}
