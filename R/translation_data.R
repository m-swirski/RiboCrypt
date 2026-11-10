# Translation stores are explicit, versioned inputs; no implicit data discovery.
translation_option <- function(browser_options, key, default = '') {
  # Missing [[ keys throw for named character vectors, the standard app format.
  if (is.null(browser_options) || !key %in% names(browser_options)) return(default)
  value <- browser_options[[key]]
  if (length(value) != 1L || is.na(value) || !nzchar(value)) return(default)
  as.character(value)
}

translation_data_dir <- function(browser_options = NULL) {
  path.expand(translation_option(browser_options, 'translation_data_dir',
                                Sys.getenv('RIBOCRYPT_TRANSLATION_DIR', '')))
}

translation_read_bundle <- function(directory) {
  path <- file.path(directory, 'atlas.rds')
  if (!file.exists(path)) stop('Configure translation_data_dir with a prepared Translation dataset.', call. = FALSE)
  x <- readRDS(path)
  if (!identical(x$version, 1L)) stop('Unsupported Translation data version.')
  stopifnot(identical(rownames(x$TE), rownames(x$RNA)),
            identical(colnames(x$TE), colnames(x$RNA)), !anyDuplicated(x$genes$gene_id),
            identical(rownames(x$TE), x$genes$gene_id), all(is.finite(x$TE)), all(is.finite(x$RNA)))
  x
}

translation_neighbours <- function(bundle, gene_id, contexts = colnames(bundle$TE), n = 50L, precomputed = NULL) {
  if (!gene_id %in% rownames(bundle$TE)) stop('Gene is not in this atlas.')
  if (length(contexts) < 3L || anyDuplicated(contexts) || !all(contexts %in% colnames(bundle$TE)))
    stop('Select at least three distinct supported contexts for an exploratory neighbour query.')
  score <- function(x) {
    x <- x[, contexts, drop = FALSE]; z <- x - rowMeans(x); g <- match(gene_id, rownames(z))
    v <- rowSums(z^2); den <- v + v[g]
    r <- 2 * drop(z %*% z[g, ]) / den
    r[!is.finite(r) | v <= 1e-8 | v[g] <= 1e-8] <- NA_real_; r[g] <- NA_real_
    r
  }
  rna <- score(bundle$RNA)
  if (!is.null(precomputed) && identical(contexts,colnames(bundle$TE)) && n<=50L) {
    ix_pre<-which(precomputed$query==gene_id);z<-precomputed[ix_pre,];z<-z[order(abs(z$TE_rho),decreasing=TRUE),]
    z<-head(z,n);ix<-match(z$gene_id,bundle$genes$gene_id)
    if(anyNA(ix))stop('Precomputed neighbours and gene universe differ.')
    te<-rep(NA_real_,nrow(bundle$TE));te[ix]<-z$TE_rho
  } else {
    te <- score(bundle$TE);ix <- head(order(abs(te), decreasing = TRUE, na.last = NA), n)
  }
  data.table::data.table(gene_id=bundle$genes$gene_id[ix], symbol=bundle$genes$symbol[ix],
                        TE_rho=te[ix], RNA_rho=rna[ix], contexts=length(contexts),
                        direction=ifelse(te[ix] > 0, 'Coordinated', 'Opposing'))
}

# Convert Run membership into assay-library-level biological units. Complete
# technical sets are mandatory, and independent library preparations are not pooled.
translation_pair_store <- function(store, mapping) {
  mapping <- as.data.frame(mapping, stringsAsFactors=FALSE)
  if (!all(c('Run','biosample_id') %in% names(mapping))) stop('Mapping needs Run and biosample_id columns.')
  mapping$Run <- as.character(mapping$Run); mapping$biosample_id <- trimws(as.character(mapping$biosample_id))
  mapping$biosample_id[mapping$biosample_id %in% c('', 'NA', 'N/A')] <- NA_character_
  if (anyNA(mapping$Run) || any(!nzchar(mapping$Run)) || anyDuplicated(mapping$Run)) stop('Run identifiers must be nonempty and unique.')
  members <- store$run_members
  if (anyDuplicated(members$Run)) stop('Count store has ambiguous Run membership.')
  if (any(!mapping$Run %in% members$Run)) stop('Mapping contains Runs absent from the count store: ', paste(head(setdiff(mapping$Run,members$Run)),collapse=', '))
  members$biosample_id <- mapping$biosample_id[match(members$Run,mapping$Run)]
  libs <- unique(members$library_id); rows <- list()
  for (lib in libs) {
    m <- members[members$library_id == lib, , drop=FALSE]; ids <- unique(m$biosample_id[!is.na(m$biosample_id)])
    if (!length(ids)) next
    if (length(ids) != 1L || anyNA(m$biosample_id)) stop('Incomplete or conflicting technical Run mapping for ', lib)
    a <- unique(m$assay); if (length(a) != 1L || !a %in% c('RNA','RIBO')) stop('Invalid assay library.')
    rows[[length(rows)+1L]] <- data.frame(biosample_id=ids,assay=a,library_id=lib,technical_runs=nrow(m))
  }
  if (!length(rows)) stop('No mapped assay libraries.')
  inventory <- do.call(rbind,rows)
  if (anyDuplicated(paste(inventory$biosample_id,inventory$assay)))
    stop('Multiple independent assay libraries map to one biosample. Review library preparation/replicate identities before pooling.')
  units <- sort(intersect(inventory$biosample_id[inventory$assay=='RNA'],inventory$biosample_id[inventory$assay=='RIBO']))
  if (!length(units)) stop('No biosamples have both RNA and Ribo libraries.')
  select <- function(assay, counts) {
    r <- inventory[inventory$assay==assay,]; ids <- r$library_id[match(units,r$biosample_id)]
    if (!all(ids %in% colnames(counts))) stop('Mapped libraries lack count columns.')
    x <- counts[,ids,drop=FALSE];colnames(x)<-units;x
  }
  R <- select('RIBO',store$Ribo); A <- select('RNA',store$RNA)
  if (!identical(rownames(R),rownames(A))) stop('RNA/Ribo feature ordering differs.')
  md <- store$units
  # Imported IDs can differ from defaults: recover metadata through library identity.
  rlib <- inventory$library_id[match(paste(units,'RIBO'),paste(inventory$biosample_id,inventory$assay))]
  meta <- md[match(rlib,md$ribo_library_id),,drop=FALSE]
  meta$biosample_id <- units
  meta$rna_library_id <- inventory$library_id[match(paste(units,'RNA'),paste(inventory$biosample_id,inventory$assay))]
  for (col in intersect(c('condition','study','replicate','match_level'),names(mapping))) {
    vals <- vapply(units,function(id) {
      z <- unique(as.character(mapping[mapping$biosample_id %in% id,col])); z <- z[!is.na(z)&nzchar(z)]
      if (length(z)>1L) stop('Conflicting ',col,' for biosample ',id)
      if(length(z)) z else NA_character_
    },character(1))
    use <- !is.na(vals);meta[[col]][use]<-vals[use]
  }
  list(Ribo=R,RNA=A,units=meta,inventory=inventory,
       unmatched_biosamples=setdiff(unique(inventory$biosample_id),units),run_members=members,genes=store$genes,
       count_region=store$count_region,reference=store$reference)
}

#' Build Translation counts from Run-resolved ORFik collections
#' @param ribo_exp,rna_exp ORFik experiments preserving individual Runs.
#' @param run_map Reviewed table containing Run and biosample_id, and optionally
#'   condition, study, replicate and match_level.
#' @param ribo_counts,rna_counts Compatible unnormalized feature-by-Run matrices.
#'   Defaults import regional collection count tables. Supply compatible feature
#'   matrices if the collection annotations/count regions differ.
#' @return Validated paired count matrices and biological-unit metadata.
#' @export
ribocrypt_translation_pairs <- function(ribo_exp, rna_exp, run_map,
    ribo_counts=NULL, rna_counts=NULL) {
  extract <- function(df,counts,region) {
    runs <- as.character(df$Run)
    if (!length(runs) || anyNA(runs) || anyDuplicated(runs)) stop('Collection must preserve distinct Run identifiers.')
    if (is.null(counts)) {
      z <- ORFik::countTable(df,region,type='summarized')
      counts <- SummarizedExperiment::assay(z)
      cd <- as.data.frame(SummarizedExperiment::colData(z))
      if (!'Run' %in% names(cd)) stop('Collection count colData must retain Run; provide explicitly Run-keyed matrices otherwise.')
      if (anyDuplicated(cd$Run) || !all(runs %in% cd$Run)) stop('Collection count Run identities are incomplete or ambiguous.')
      counts <- counts[,match(runs,cd$Run),drop=FALSE];colnames(counts)<-runs
    }
    if (is.null(colnames(counts)) || !all(runs %in% colnames(counts))) stop('Counts must be keyed by Run.')
    counts[,runs,drop=FALSE]
  }
  if (!identical(BiocGenerics::organism(ribo_exp), BiocGenerics::organism(rna_exp))) stop('Collection species differ.')
  if (!identical(ORFik::findFa(ribo_exp),ORFik::findFa(rna_exp))) stop('Collection genome references differ.')
  if (methods::is(ribo_exp,'experiment') && methods::is(rna_exp,'experiment') && !identical(ribo_exp@txdb,rna_exp@txdb))stop('Collection transcript annotations differ.')
  R<-extract(ribo_exp,ribo_counts,'cds');A<-extract(rna_exp,rna_counts,'mrna')
  if (is.null(rownames(R)) || is.null(rownames(A)) || anyDuplicated(rownames(R)) || anyDuplicated(rownames(A))) stop('Provide unique compatible feature identities in both assays.')
  common<-intersect(rownames(R),rownames(A));if(!length(common))stop('Count regions have no compatible feature identities.')
  R<-R[match(common,rownames(R)),,drop=FALSE];A<-A[match(common,rownames(A)),,drop=FALSE]
  rm<-data.frame(Run=colnames(R),library_id=colnames(R),assay='RIBO')
  am<-data.frame(Run=colnames(A),library_id=colnames(A),assay='RNA')
  units<-data.frame(ribo_library_id=colnames(R),biosample_id=colnames(R),condition=as.character(ribo_exp$condition),study=ORFik::name(ribo_exp),replicate=as.character(ribo_exp$rep),match_level='user_review_required')
  store<-list(Ribo=R,RNA=A,run_members=rbind(rm,am),units=units,
              genes=data.frame(gene_id=rownames(R),symbol=rownames(R)),count_region='Ribo CDS / RNA mRNA',reference=ORFik::findFa(ribo_exp))
  # Collection Run columns are summed within a reviewed assay_library_id when
  # supplied; otherwise each Run is a distinct library (no assumed replication).
  if ('assay_library_id' %in% names(run_map)) {
    ids<-as.character(run_map$assay_library_id[match(store$run_members$Run,run_map$Run)])
    missing<-is.na(ids)|!nzchar(ids);ids[missing]<-store$run_members$library_id[missing]
    keys<-paste(store$run_members$assay,ids,sep='::');store$run_members$library_id<-keys
    collapse<-function(x,assay) {
      k<-keys[store$run_members$assay==assay];t(rowsum(t(x),group=k,reorder=FALSE))
    }
    store$Ribo<-collapse(R,'RIBO');store$RNA<-collapse(A,'RNA')
    old<-store$units;rkeys<-keys[seq_len(nrow(rm))]
    for (key in unique(rkeys)) {
      ix<-which(rkeys==key)
      if (length(unique(old$condition[ix]))>1L) stop('Technical library has conflicting conditions.')
    }
    store$units<-old[match(colnames(store$Ribo),rkeys),,drop=FALSE]
    store$units$ribo_library_id<-colnames(store$Ribo)
  }
  result<-translation_pair_store(store,run_map)
  result$species<-BiocGenerics::organism(ribo_exp)
  result
}

# Pair QC is distinct from biological correspondence and differential significance.
translation_pair_r2 <- function(pairs, qc = NULL) {
  R<-pairs$Ribo;A<-pairs$RNA
  if(!identical(dimnames(R),dimnames(A)))stop('Pair QC requires aligned features and units.')
  if(any(!is.finite(R))||any(!is.finite(A))||any(R<0)||any(A<0))stop('Invalid counts for pair QC.')
  frozen<-!is.null(qc)&&identical(pairs$reference,qc$reference)&&all(c(qc$feature_ids,qc$other_ids)%in%rownames(R))
  components<-function(x) {
    if(frozen)x<-rbind(x[qc$feature_ids,,drop=FALSE],others=colSums(x[qc$other_ids,,drop=FALSE]))
    z<-log(pmax(x,1));z-colMeans(z)[col(z)]
  }
  r<-components(R);a<-components(A)
  den<-sqrt(colSums(r*r)*colSums(a*a));rho<-colSums(r*a)/den
  rho[!is.finite(rho)]<-NA_real_
  list(R2=pmin(1,pmax(0,rho^2)),Pearson_r=rho,
       definition=if(frozen)qc$definition else 'Squared Pearson correlation across supplied common feature CLR-count components; zeros floored to one (reference cutoff is exploratory)')
}

translation_attach_pair_qc <- function(pairs, qc = NULL) {
  u<-pairs$units
  if(!is.null(qc)) {
    z<-qc$units[match(u$ribo_library_id,qc$units$ribo_library_id),,drop=FALSE]
    unchanged<-!is.na(z$qc_rna_library_id)&u$rna_library_id==z$qc_rna_library_id
    u$matched_R2<-ifelse(unchanged,z$matched_R2,NA_real_)
    u$matching_score<-ifelse(unchanged,z$matching_score,NA_real_)
    u$Ribo_CDS_fraction_proxy<-z$Ribo_CDS_fraction_proxy
  } else if(!'matched_R2'%in%names(u))u$matched_R2<-NA_real_
  missing<-!is.finite(u$matched_R2)
  if(any(missing)) {
    sub<-pairs;sub$Ribo<-pairs$Ribo[,missing,drop=FALSE];sub$RNA<-pairs$RNA[,missing,drop=FALSE]
    score<-translation_pair_r2(sub,qc);u$matched_R2[missing]<-score$R2
    pairs$pair_qc_definition<-score$definition
  } else pairs$pair_qc_definition<-qc$definition
  pairs$units<-u;pairs
}

translation_filter_pair_qc <- function(units, minimum = .188, enabled = TRUE, include_unscored = FALSE) {
  if(!isTRUE(enabled))return(units)
  if(length(minimum)!=1L||!is.finite(minimum)||minimum<0||minimum>1)stop('Pair R² cutoff must be between zero and one.')
  score<-if('matched_R2'%in%names(units))units$matched_R2 else rep(NA_real_,nrow(units))
  keep<-is.finite(score)&score>=minimum
  if(isTRUE(include_unscored))keep<-keep|!is.finite(score)
  units[keep,,drop=FALSE]
}
