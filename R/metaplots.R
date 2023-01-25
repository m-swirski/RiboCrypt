unlistToExtremities <- function(grl) {
  gr <- unlistGrl(grl)
  gr$names <- names(gr)
  dt <- as.data.table(gr)
  dt <- dt[,.(seqnames = seqnames[1], start = min(start), end = max(end), strand = strand[1]), by = names]
  return(GRanges(dt))
}

#' Distance to following range
#'
#' @param grl a GRangesList
#' @param grl2 a GRangesList, default 'grl'
#' @param ignore.strand logical, default FALSE
#' @importFrom IRanges precede follow distance
#' @return numeric vector of distance
distanceToFollowing <- function(grl, grl2 = grl, ignore.strand = FALSE) {
  stops <- stopRegion(grl,downstream = 0, upstream = 0) %>% unlistGrl
  starts <- grl2 %>% unlistToExtremities()
  if (!ignore.strand){
    followers <- precede(stops,starts)
  } else {
    minus_strand <- which(stops@strand == "-")
    plus_strand <- which(stops@strand == "+")
    followers_plus <- precede(stops[plus_strand],starts, ignore.strand = TRUE)
    followers_minus <- follow(stops[minus_strand],starts, ignore.strand = TRUE)
    followers <- c()
    if (length(minus_strand > 0)) followers[minus_strand] <- followers_minus
    if (length(plus_strand > 0)) followers[plus_strand] <- followers_plus
  }
  nas <- is.na(followers)
  dists <- distance(stops[!is.na(followers)], starts[followers[!is.na(followers)]], ignore.strand = ignore.strand)
  out <- 1:length(grl)
  out[!nas] <- dists
  out[nas] <- NA
  out
}

distanceToPreceding <- function(grl, grl2 = grl, ignore.strand = FALSE) {
  stops <- grl2 %>% unlistToExtremities()
  starts <- startRegion(grl,downstream = 0,upstream = 0) %>% unlistGrl()
  if (!ignore.strand){
    preceders <- follow(starts,stops)
  } else {
    minus_strand <- which(starts@strand == "-")
    plus_strand <- which(starts@strand == "+")
    preceders_plus <- follow(starts[plus_strand],stops, ignore.strand = TRUE)
    preceders_minus <- precede(starts[minus_strand],stops, ignore.strand=TRUE)
    preceders <- c()
    if (length(minus_strand > 0)) preceders[minus_strand] <- preceders_minus
    if (length(plus_strand > 0)) preceders[plus_strand] <- preceders_plus
  }
  nas <- is.na(preceders)
  dists <- distance(starts[!is.na(preceders)], stops[preceders[!is.na(preceders)]], ignore.strand = ignore.strand)
  out <- 1:length(grl)
  out[!nas] <- dists
  out[nas] <- NA
  out

}

extendTrailersUntil <- function(grl, grl2=grl, extension = 500, until = 200, min_ext = 25, ...) {
  dists <- distanceToFollowing(grl,grl2, ...)
  dists[is.na(dists)] <- extension + until + 1
  diff <- pmax(dists - until, min_ext)
  extension <- pmin(extension, diff)
  grl <- extendTrailers(grl, extension)
  return(grl)
}

extendLeadersUntil <- function(grl, grl2=grl, extension = 500, until = 200, min_ext = 25, ...) {
  dists <- distanceToPreceding(grl,grl2, ...)
  dists[is.na(dists)] <- extension + until + 1
  diff <- pmax(dists - until,min_ext)
  extension <- pmin(extension, diff)


  grl <- extendLeaders(grl, extension)
  return(grl)
}

getStopWindow <- function(grl, upstream = 200, downstream = 500, min_upstream_dist = 50, min_downstream_dist = 200, min_ext = 25,...) {
  ext_grl <- grl %>% stopRegion(upstream = 0, downstream = 0) %>% extendTrailersUntil(.,grl,extension=downstream,until=min_downstream_dist,min_ext=min_ext,...)

  ext_grl <- extendLeadersUntil(ext_grl, startRegion(grl), extension = upstream, until = min_upstream_dist, min_ext = min_ext,...)
  return(ext_grl)

}

getStartWindow <- function(grl, upstream = 500, downstream = 200, min_upstream_dist = 200, min_downstream_dist = 50, min_ext = 25,...) {
  ext_grl <- startRegion(grl, upstream = 0, downstream = 0) %>% extendLeadersUntil(.,grl,extension=upstream,until=min_upstream_dist,min_ext=min_ext,...)

  ext_grl <- extendTrailersUntil(ext_grl, stopRegion(grl, upstream = 0, downstream = 0), extension = downstream, until = min_downstream_dist, min_ext = min_ext,...)
  return(ext_grl)
}

stopCoverage <- function(reads, grl, upstream = 200, downstream = 500,min_upstream_dist = 50, min_downstream_dist = 200, min_ext = 25,...) {
  stopWindow <- getStopWindow(grl, upstream,downstream,min_upstream_dist,min_downstream_dist,min_ext,...)
  numbering <- distance(grl %>% stopRegion(upstream = 0, downstream = 0) %>% unlistGrl, startRegion(stopWindow,upstream = 0, downstream = 0) %>% unlistGrl)
  cov <- coveragePerTiling(stopWindow, reads, as.data.table = TRUE)
  numbering <- numbering + 2
  genes <- NULL # avoid BiocCheck warning
  cov[,position := position - numbering[genes]]
  return(cov)
}
startCoverage <- function(reads, grl, upstream = 500, downstream = 200, min_upstream_dist = 200, min_downstream_dist = 50, min_ext = 25,...) {
  startWindow <- getStartWindow(grl, upstream,downstream,min_upstream_dist,min_downstream_dist,min_ext,...)
  numbering <- distance(startRegion(grl, upstream = 0, downstream = 0) %>% unlistGrl, startRegion(startWindow,upstream = 0, downstream = 0) %>% unlistGrl)
  cov <- coveragePerTiling(startWindow, reads, as.data.table = TRUE)
  numbering <- numbering + 2
  genes <- NULL # avoid BiocCheck warning
  cov[,position := position - numbering[genes]]
  return(cov)
}

getMetaCoverage <- function(reads, grl, outward = 500, inward = 200, min_outward_distance = 200, min_inward_distance = 50, min_ext = 25, transcriptNormalize = TRUE, withFrames = TRUE,...) {
  start_cov <- startCoverage(reads, grl, upstream = outward, downstream = inward, min_upstream_dist = min_outward_distance, min_downstream_dist = min_inward_distance, min_ext = min_ext,...)
  stop_cov <- stopCoverage(reads,grl, upstream = inward, downstream = outward , min_upstream_dist = min_inward_distance, min_downstream_dist = min_outward_distance, min_ext = min_ext,...)
  genes <- type <- count <- NULL # avoid BiocCheck warning
  if (transcriptNormalize) {
    # total_counts <- rbind(start_cov,stop_cov)[,.(count = sum(count)), by = genes]$count
    start_wind <- getStartWindow(grl, upstream = outward + as.integer(outward/2), downstream = inward, min_upstream_dist = min_outward_distance, min_downstream_dist = min_inward_distance, min_ext = min_ext,...)
    stop_wind <- getStopWindow(grl, upstream = inward, downstream = outward + as.integer(outward/2), min_upstream_dist = min_inward_distance, min_downstream_dist = min_outward_distance, min_ext = min_ext,...)
    start_counts <- countOverlapsW(start_wind, reads, weight = "score")
    stop_counts <- countOverlapsW(stop_wind, reads, weight = "score")
    total_counts <- start_counts + stop_counts
    total_counts <- pmax(total_counts, 1)
    start_cov[, count := count / total_counts[genes]]
    stop_cov[, count := count / total_counts[genes]]
  }

  start_meta <- start_cov[,.(count = sum(count)), by = position]
  stop_meta <- stop_cov[,.(count = sum(count)), by = position]

  start_meta <- start_meta[order(position)]
  stop_meta <- stop_meta[order(position)]

  start_meta[,type := "start"]
  stop_meta[,type := "stop"]
  all_meta <- rbind(start_meta, stop_meta)
  all_meta[, index := 1:.N]
  if (withFrames) all_meta[,frames := as.factor(position %% 3)]
  return(all_meta)

}

metaPlot <- function(reads, grl, outward = 500, inward = 200, min_outward_distance = 200, min_inward_distance = 50, min_ext = 25, transcriptNormalize = TRUE, withFrames = TRUE, col = "black",...) {
  grl <- grl %>% extendTrailers(1)
  metaCov <- getMetaCoverage(reads, grl, outward, inward, min_outward_distance, min_inward_distance, min_ext, transcriptNormalize=transcriptNormalize, withFrames=withFrames,...)
  mybreaks <- c(outward/2 + 1, outward + 1, outward + inward/2 + 1, outward + 1.5 * inward + 2,outward + 2*inward + 2, 1.5 * outward + 2 * inward + 2)
  myticks <- metaCov$position[mybreaks]
  count <- NULL # Avoid warning
  if (withFrames) {
    ggplot(metaCov) +
      geom_vline(aes(xintercept = inward + outward)) + 
      geom_col(aes(y = count, x = index, color = frames), size = 0.75) +
      theme(legend.position = "none") +
      ylab("transcript normalized coverage") +
      xlab("relative position [nt]") +
      theme(plot.margin = unit(c(0,0,0,0), "pt")) +
      scale_x_continuous(expand = c(0,0), breaks = mybreaks, labels = myticks) +
      theme_bw()
  } else {
    ggplot(metaCov)  +
      geom_area(aes(y = count, x = index), fill = col, position = "identity") +
      theme(plot.margin = unit(c(0,0,0,0), "pt"))+
      scale_x_continuous(expand = c(0,0), breaks = mybreaks, labels = myticks)
  }
}

threePlots <- function(TSS,PAS,RNA, grl, outward = 500, inward = 200, min_outward_distance = 200, min_inward_distance = 50, min_ext = 25, transcriptNormalize = TRUE, withFrames = TRUE, col = "black",...) {
  grl <- grl %>% extendTrailers(1)
  metaCovTSS <- getMetaCoverage(TSS, grl, outward, inward, min_outward_distance, min_inward_distance, min_ext, transcriptNormalize=transcriptNormalize, withFrames=FALSE)[,seqType := 'TSS']
  metaCovPAS <- getMetaCoverage(PAS, grl, outward, inward, min_outward_distance, min_inward_distance, min_ext, transcriptNormalize=transcriptNormalize, withFrames=FALSE)[,seqType := 'PAS']
  metaCovRNA <- getMetaCoverage(RNA, grl, outward, inward, min_outward_distance = 0, min_inward_distance = 0 , min_ext = 0, transcriptNormalize=transcriptNormalize, withFrames=FALSE,ignore.strand = TRUE)[,seqType := 'RNA']
  metaCov <- rbind(metaCovTSS,metaCovPAS, metaCovRNA)
  mybreaks <- c(outward/2 + 1, outward + 1, outward + inward/2 + 1, outward + 1.5 * inward + 2,outward + 2*inward + 2, 1.5 * outward + 2 * inward + 2)
  myticks <- metaCovTSS$position[mybreaks]
  endlines <- metaCovTSS[position == 0]$index
  midline <- nrow(metaCovTSS) / 2
  count <- NULL # Avoid warning
  metaCov[,count := count/max(count), by = seqType]



  ggplot()  +
    geom_area(data = metaCov[seqType == 'RNA'],mapping = aes(y = count, x = index), fill = "grey35", alpha = 0.8, position = "identity") +
    geom_area(data = metaCov[seqType == 'TSS'], mapping = aes(y = count, x = index), fill = "orange", alpha = 0.6, position = "identity") +
    geom_area(data = metaCov[seqType == 'PAS'],mapping = aes(y = count, x = index), fill = "purple", alpha = 0.6, position = "identity") +
    geom_vline(xintercept = endlines, linetype = 2) +
    geom_vline(xintercept = midline) +
    theme(plot.margin = unit(c(0,0,0,0), "pt"))+
    scale_x_continuous(expand = c(0,0), breaks = mybreaks, labels = myticks)

}

fourPlots <- function(ribo, TSS,PAS,RNA, grl, outward = 500, inward = 200, min_outward_distance = 200, min_inward_distance = 50, min_ext = 25, transcriptNormalize = TRUE, withFrames = TRUE, col = "black",...) {
  grl <- grl %>% extendTrailers(1)
  metaCovRibo <-   getMetaCoverage(ribo, grl, outward, inward, min_outward_distance, min_inward_distance, min_ext, transcriptNormalize=transcriptNormalize, withFrames=withFrames)[,seqType := 'RIBO']
  metaCovTSS <- getMetaCoverage(TSS, grl, outward, inward, min_outward_distance, min_inward_distance, min_ext, transcriptNormalize=transcriptNormalize, withFrames=withFrames)[,seqType := 'TSS']
  metaCovPAS <- getMetaCoverage(PAS, grl, outward, inward, min_outward_distance, min_inward_distance, min_ext, transcriptNormalize=transcriptNormalize, withFrames=withFrames)[,seqType := 'PAS']
  metaCovRNA <- getMetaCoverage(RNA, grl, outward, inward, min_outward_distance=0, min_inward_distance=0, min_ext=0, transcriptNormalize=transcriptNormalize, withFrames=withFrames,...)[,seqType := 'RNA']
  metaCov <- rbind(metaCovRibo, metaCovTSS,metaCovPAS, metaCovRNA)
  mybreaks <- c(outward/2 + 1, outward + 1, outward + inward/2 + 1, outward + 1.5 * inward + 2,outward + 2*inward + 2, 1.5 * outward + 2 * inward + 2)
  count <- NULL # Avoid warning
  myticks <- metaCovTSS$position[mybreaks]
  endlines <- metaCovTSS[position == 0]$index
  midline <- nrow(metaCovTSS) / 2
  metaCov[,count := count/max(count), by = seqType]

  rna_seqs <-  ggplot()  +
    geom_area(data = metaCov[seqType == 'RNA'],mapping = aes(y = count, x = index), fill = "grey35", alpha = 0.8, position = "identity") +
    geom_area(data = metaCov[seqType == 'TSS'], mapping = aes(y = count, x = index), fill = "orange", alpha = 0.6, position = "identity") +
    geom_area(data = metaCov[seqType == 'PAS'],mapping = aes(y = count, x = index), fill = "purple", alpha = 0.6, position = "identity") +
    geom_vline(xintercept = endlines, linetype = 2) +
    geom_vline(xintercept = midline) +
    theme(plot.margin = unit(c(0,0,0,0), "pt"))+
    scale_x_continuous(expand = c(0,0), breaks = mybreaks, labels = myticks) +
    theme(plot.margin = unit(c(0,0,0,0), "pt")) +
    scale_x_continuous(expand = c(0,0)) +
    theme(axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank()) +
    theme(legend.position = "none") +
    ylab(NULL)

  ribo_seq <- ggplot()  +
    geom_vline(xintercept = endlines, linetype = 2) +
    geom_vline(xintercept = midline) +
    geom_line(data = metaCov[seqType == 'RIBO'], mapping = aes(y = count, x = index, color = frames), size = 0.75) +
    theme(plot.margin = unit(c(0,0,0,0), "pt"))+
    scale_x_continuous(expand = c(0,0), breaks = mybreaks, labels = myticks) +
    theme(legend.position = "none") +
    ylab(NULL)


  multiomics_plot <- subplot(rna_seqs %>% automateTicksRNA(),
                             ribo_seq %>% automateTicks(),
                             margin = 0,
                             nrows = 2,
                             heights = c(0.5,0.5),
                             shareX = TRUE,
                             titleY = TRUE,
                             titleX = TRUE)
  multiomics_plot <- multiomics_plot %>% plotly::config(
    toImageButtonOptions = list(format = "svg"))

  return(multiomics_plot)

}

fivePlots <- function(ribo, TSS,PAS,RNA, grl, outward = 500, inward = 200, min_outward_distance = 200, min_inward_distance = 50, min_ext = 25, transcriptNormalize = TRUE, withFrames = TRUE, col = "black",...) {
  grl <- grl %>% extendTrailers(1)
  metaCovRibo <-   getMetaCoverage(ribo, grl, outward, inward, min_outward_distance, min_inward_distance, min_ext, transcriptNormalize=transcriptNormalize, withFrames=withFrames)[,seqType := 'RIBO']
  metaCovTSS <- getMetaCoverage(TSS, grl, outward, inward, min_outward_distance, min_inward_distance, min_ext, transcriptNormalize=transcriptNormalize, withFrames=withFrames)[,seqType := 'TSS']
  metaCovPAS <- getMetaCoverage(PAS, grl, outward, inward, min_outward_distance, min_inward_distance, min_ext, transcriptNormalize=transcriptNormalize, withFrames=withFrames)[,seqType := 'PAS']
  metaCovRNA <- getMetaCoverage(RNA, grl, outward, inward, min_outward_distance=0, min_inward_distance=0, min_ext=0, transcriptNormalize=transcriptNormalize, withFrames=withFrames,...)[,seqType := 'RNA']

  metaCov <- rbind(metaCovRibo, metaCovTSS,metaCovPAS, metaCovRNA)


  metaCovReadlength <- readLengthMeta(ribo,grl, outward = outward, inward = inward, min_outward_distance = min_outward_distance, min_inward_distance = min_inward_distance, min_ext = min_ext, transcriptNormalize = transcriptNormalize)

  metaCovReadlength$index <- 1:nrow(metaCovReadlength)
  count <- NULL # Avoid warning
  mybreaks <- c(outward/2 + 1, outward + 1, outward + inward/2 + 1, outward + 1.5 * inward + 2,outward + 2*inward + 2, 1.5 * outward + 2 * inward + 2)
  myticks <- metaCovTSS$position[mybreaks]
  endlines <- metaCovTSS[position == 0]$index
  midline <- nrow(metaCovTSS) / 2
  metaCov[,count := count/max(count), by = seqType]


  rna_seqs <-  ggplot()  +
    geom_area(data = metaCov[seqType == 'RNA'],mapping = aes(y = count, x = index), fill = "grey35", alpha = 0.8, position = "identity") +
    geom_area(data = metaCov[seqType == 'TSS'], mapping = aes(y = count, x = index), fill = "orange", alpha = 0.6, position = "identity") +
    geom_area(data = metaCov[seqType == 'PAS'],mapping = aes(y = count, x = index), fill = "purple", alpha = 0.6, position = "identity") +
    geom_vline(xintercept = endlines, linetype = 2) +
    geom_vline(xintercept = midline) +
    theme(plot.margin = unit(c(0,0,0,0), "pt"))+
    scale_x_continuous(expand = c(0,0), breaks = mybreaks, labels = myticks) +
    theme(plot.margin = unit(c(0,0,0,0), "pt")) +
    scale_x_continuous(expand = c(0,0)) +
    theme(axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank()) +
    theme(legend.position = "none") +
    ylab(NULL)

  ribo_seq <- ggplot()  +
    geom_vline(xintercept = endlines, linetype = 2) +
    geom_vline(xintercept = midline) +
    geom_line(data = metaCov[seqType == 'RIBO'], mapping = aes(y = count, x = index, color = frames), size = 0.75) +
    theme(plot.margin = unit(c(0,0,0,0), "pt"))+
    scale_x_continuous(expand = c(0,0), breaks = mybreaks, labels = myticks) +
    theme(legend.position = "none") +
    ylab(NULL)


  ratio_plot <- ggplot(metaCovReadlength) +
    geom_vline(xintercept = endlines, linetype = 2) +
    geom_vline(xintercept = midline) +
    geom_line(aes(x = index, y =ratio * count), color = "magenta", size = 0.75) +
    theme(plot.margin = unit(c(0,0,0,0), "pt"))+
    scale_x_continuous(expand = c(0,0), breaks = mybreaks, labels = myticks) +
    theme(legend.position = "none") +
    ylab(NULL)


  multiomics_plot <- subplot(rna_seqs %>% ggplotly,
                             ribo_seq %>% ggplotly,
                             ratio_plot %>% ggplotly,
                             margin = 0,
                             nrows = 3,
                             heights = c(0.33,0.33,0.33),
                             shareX = TRUE,
                             titleY = TRUE,
                             titleX = TRUE)
  multiomics_plot <- multiomics_plot %>% plotly::config(
    toImageButtonOptions = list(format = "svg"))

  return(multiomics_plot)
}
