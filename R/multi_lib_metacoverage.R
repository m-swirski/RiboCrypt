multiLib_to_metaWindow <- function(collection_df,region = c("stop","start")[1], windowUpstream = 30, windowDownstream = 30, minFiveUTR = 30L, minCDS = 150L, minThreeUTR = 30L, cores = 1, scoring = "zscore"){
  fps <- sapply(1:nrow(collection_df),function(x) try(filepath(collection_df[x,], type = "cov")))
  cov_exist <- fps != "Error in FUN(X[[i]], ...) : \n  File did not exist, did you create covRle yet?\n" 
  fps <- fps[cov_exist]
  
  minFiveUTR <- max(minFiveUTR, windowUpstream + 1)
  minThreeUTR <- max(minThreeUTR, windowDownstream + 1)
  trxs <- filterTranscripts(collection_df,minFiveUTR = minFiveUTR, minCDS = minCDS, minThreeUTR = minThreeUTR, longestPerGene = TRUE)
  
  cds <- loadRegion(collection_df,"cds", names.keep = trxs)
  trx <- loadRegion(collection_df, "transcript", names.keep = trxs)
  browser()
  if (region == "stop") {
    windows <- stopRegion(cds, tx = trx, upstream = windowUpstream, downstream = windowDownstream)
  } else {
    windows <- startRegion(cds, tx = trx, upstream = windowUpstream, downstream = windowDownstream)
  }
  output <- mclapply(fps, function(x) coverageScorings(coveragePerTiling(windows, fimport(x), as.data.table = TRUE)[,position := -windowDownstream:windowUpstream, by = genes], scoring = scoring, copy.dt = FALSE ) , mc.cores= cores)
  names(output) <- bamVarName(collection_df)[cov_exist]
  output <- rbindlistst(output, idcol = "library")
  return(output)
}

collection_to_metawindow <- function(collection_df, region = c("stop","start")[1], windowUpstream = 30, windowDownstream = 30, minFiveUTR = 30L, minCDS = 150L, minThreeUTR = 30L, scoring = "zscore", cores = 1){
  collection_path <- file.path(resFolder(collection_df), "collection_tables")
  minFiveUTR <- max(minFiveUTR, windowUpstream + 1)
  minThreeUTR <- max(minThreeUTR, windowDownstream + 1)
  trxs <- filterTranscripts(collection_df,minFiveUTR = minFiveUTR, minCDS = minCDS, minThreeUTR = minThreeUTR, longestPerGene = TRUE)
  true_collection <- system(paste0("ls ", collection_path), intern = TRUE ) %>% sub(".fst","",.)
  trxs <- trxs[trxs %in% true_collection]
  
  
  cds <- loadRegion(collection_df,"cds", names.keep = trxs)
  trx <- loadRegion(collection_df, "transcript", names.keep = trxs)
  if (region == "stop") {
    windows <- stopRegion(cds, tx = trx, upstream = windowUpstream, downstream = windowDownstream)
  } else {
    windows <- startRegion(cds, tx = trx, upstream = windowUpstream, downstream = windowDownstream)
  }
  relcor <- pmapToTranscriptF(windows, trx) 
  relstart <- unlist(start(relcor))
  relend <- unlist(end(relcor))
  
  output <- mcmapply(function(x,start,end)  internalGetterSubsetter(x,start,end, collection_path = collection_path, windowUpstream = windowUpstream, windowDownstream = windowDownstream),
                     trxs, relstart, relend, mc.cores = cores, SIMPLIFY = FALSE)
  output <- rbindlist(output, id = "genes")
  
  output <- output[ ,coverageScorings(.SD, scoring = scoring), by = library]
  return(output)
}
internalGetterSubsetter <- function(x,start,end, collection_path, windowUpstream, windowDownstream) {
  ff <- as.data.table(fst::read.fst(paste0(collection_path,"/", x,".fst")))
  ff <- ff[, .SD[start:end,] , by = library][, position := -windowUpstream:windowDownstream,by=library ]
  return(ff)
}
