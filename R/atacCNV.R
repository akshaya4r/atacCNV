#' @import stats
#' @import GenomicRanges
#' @import plyranges
#' @import S4Vectors
#' @import BiocGenerics
#' @import BSgenome.Hsapiens.UCSC.hg38
#' @import data.table
#' @import Rsamtools
#' @import GenomeInfoDb
#' @import Biostrings
#' @import Matrix
#' @import ggplot2
#' @import cowplot
#' @import SummarizedExperiment
#' @import magrittr
#' @importFrom GenomicAlignments summarizeOverlaps
#' @return \code{NULL}
#' @export


# Imports:
# GenomicAlignments (>= 1.22.1),
# Matrix (>= 1.2.18),
# SummarizedExperiment (>= 1.16.1),
# data.table (>= 1.12.8),
# plyranges (>= 1.6.10),
# Rsamtools (>= 2.2.1),
# GenomeInfoDb (>= 1.22.0),
# BSgenome.Hsapiens.UCSC.hg38 (>= 1.4.1),
# GenomicRanges (>= 1.38.0),
# stats (>= 3.6.1),
# Biostrings (>= 2.54.0),
# magrittr (>= 1.5),
# BiocGenerics (>= 0.32.0),
# S4Vectors (>= 0.24.2)

atacCNV <- function(input, outdir, blacklist, windowSize, genome="BSgenome.Hsapiens.UCSC.hg38",
                    test='AD', reuse.existing=FALSE, exclude=NULL, readout="ATAC",
                    uq=0.8, lq=0.5, somyl=0.2, somyu=0.8, title_karyo=NULL, minFrags = 20000,
                    gene.annotation=NULL, threshold_blacklist_bins=0.85){

  if(reuse.existing==FALSE){
    print("Removing old file from the output folder")
    file.remove(list.files(outdir, full.names=TRUE))
  }

  if(!file.exists(file.path(outdir,"count_summary.rds"))) {
    blacklist <- read_bed(blacklist)
    windows <- makeWindows(genome = genome, blacklist = blacklist, windowSize, exclude=exclude)

    if(file_test("-d", input)){
      print("Obtaining bam file list")
      bamfiles <- Rsamtools::BamFileList(list.files(input, pattern = ".bam$", full.names = TRUE), yieldSize=100000)
      print(bamfiles)
      counts <- generateCountMatrix(bamfiles, windows)
    }
    else if(file_test("-f", input)){
      if(grepl("\\.tsv$", input)){
        print("Obtaining the fragments tsv file")
        file_fragments <- fread(input)
        colnames(file_fragments) <- c('seqnames','start','end','barcode','pcr')
        fragments <- as_granges(file_fragments)
        print(head(fragments))
        print(names(fragments))
      } else if(grepl("\\.bed$", input)){
        print("Obtaining the fragments bed file")
        fragments <- read_bed(input)
        names(mcols(fragments)) <- 'barcode'
      } else{
        stop("Please provide a fragments .tsv/.bed or a path to the directory containing all the bam files")
      }
      counts <- generateCountMatrix(fragments, windows, by="barcode", minFrags = minFrags)
    }
    saveRDS(counts, file.path(outdir,"count_summary.rds"))
  }

  counts <- readRDS(file.path(outdir,"count_summary.rds"))
  peaks <- as.data.table(assays(counts)$counts)
  colnames(peaks) <- paste0('cell-', colnames(peaks))
  rowinfo <- as.data.table(rowRanges(counts))
  peaks <- cbind(rowinfo, peaks)

  if(!file.exists(file.path(outdir,"counts_gc_corrected.rds"))) {
    corrected_counts <- peaks[, lapply(.SD, function(x) {
      fit <- stats::loess(x ~ peaks$GC)
      correction <- mean(x) / fit$fitted
      as.integer(round(x * correction))
    }), .SDcols = patterns("cell-")]
    saveRDS(corrected_counts, file.path(outdir,"counts_gc_corrected.rds"))
  }

  corrected_counts <- readRDS(file.path(outdir,"counts_gc_corrected.rds"))
  peaks <- cbind(rowinfo, corrected_counts)

  if(!is.null(gene.annotation)) {
    rowinfo.gr <- rowRanges(counts)
    rowinfo.gr <- addExpressionFactor(rowinfo.gr, gene.annotation)
    print(rowinfo.gr)
    ## corrected_counts <- corrected_counts/rowinfo.gr$numgenes
    # corrected_counts <- corrected_counts[, lapply(.SD, function(x) {
    #   # fit <- stats::loess(x ~ rowinfo.gr$numgenes)
    #   fit <- stats::loess(x ~ rowinfo.gr$genecoverage)
    #   correction <- mean(x) / fit$fitted
    #   as.integer(round(x * correction))
    # })]
    peaks <- cbind(rowinfo, corrected_counts)
    peaks <- peaks[peaks$genecoverage>0]
  }


  zeroes_per_bin <- peaks[, rowSums(.SD==0), .SDcols = patterns("cell-")]
  ncells <- length(grep("cell-", colnames(peaks)))
  peaks <- peaks[zeroes_per_bin<(threshold_blacklist_bins*ncells)]

  if(!file.exists(file.path(outdir,"results_gc_corrected.rds"))) {
    clusters_ad <- peaks[, lapply(.SD, function(x) {
      k <- 3 # k produces 2^k - 1 breakpoints
      minsize <- 5
      peaksperchrom <- split(x, peaks$seqnames)
      print("Calculating distance AD")
      results <- lapply(peaksperchrom, function(x2) {
        # x2 <- runmed(x2, k=11, endrule = 'median')
        getbp(x2, k = k, minsize = minsize, test=test)
      })
      # clusterperchrom <- lapply(results, '[[', 3)
      # cl <- 0
      # clusters <- vector()
      # for(item in clusterperchrom){
      #       item <- item + cl
      #       clusters <- append(clusters, item)
      #       cl <- cl + length(unique(item))
      # }
      # return(clusters)
    }), .SDcols = patterns("cell-")]
    saveRDS(clusters_ad, file.path(outdir, "results_gc_corrected.rds"))
  }
  print("Successfully identified breakpoints")

  names_seq <- levels(peaks$seqnames)

  clusters_ad <- readRDS(file.path(outdir,"results_gc_corrected.rds"))
  breakpoints <- lapply(clusters_ad, function(x) { lapply(x,'[[', 1) })
  distances <- lapply(clusters_ad, function(x) { lapply(x,'[[', 2) })
  result.dt <- Map(function(bp, dist){
    names(bp) <- names_seq
    names(dist) <- names_seq
    dtlist <- Map(function(per_chr_bp, per_chr_dist){
      data.table(per_chr_bp, per_chr_dist)
    }, bp, dist)
    per_cell_dt <- rbindlist(dtlist, idcol = 'Chr')
    per_cell_dt <- per_cell_dt[order(-per_chr_dist)]
  }, breakpoints, distances)

  pruned_result.dt <- lapply(result.dt, function(x){
    threshold_dist_values(x)
  })
  print("Successfully discarded irrelevant breakpoints")


  clusters_pruned <- as.data.table(Map(function(seq_data, bp){
    per_chrom_seq_data <- split(seq_data, peaks$seqnames)
    per_chrom_bp <- split(bp, bp$Chr)
    clusters_per_chrom <- Map(function(seq_data2, bp2){
      if(is.null(bp2)){
        return(rep(1, length(seq_data2)))
      } else {
        bp_to_cluster <- sort(c(1,length(seq_data2)+1,bp2$per_chr_bp))
        return(rep(1:length(diff(bp_to_cluster)),diff(bp_to_cluster)))
      }
    }, per_chrom_seq_data[names_seq], per_chrom_bp[names_seq])
    cl <- 0
    clusters <- vector()
    for(item in clusters_per_chrom){
      item <- item + cl
      clusters <- append(clusters, item)
      cl <- cl + length(unique(item))
    }
    return(clusters)
  },  peaks[, .SD, .SDcols = patterns("cell-")], pruned_result.dt))

  if(readout=="ATAC"){
    somies_ad <- Map(function(seq_data,cluster) {
      assign_gainloss(seq_data, cluster, uq=uq, lq=lq)
    }, peaks[, .SD, .SDcols = patterns("cell-")], clusters_pruned)
    print("Successfully assigned gain-loss")
    if(is.null(title_karyo)){
      title_karyo <- basename(outdir)
    }
    plot_karyo_gainloss(somies_ad = somies_ad, outdir = outdir, peaks = peaks, uq, lq, title_karyo)
    print("Successfully plotted karyogram")
  }
  if(readout=="BS"){
    somies_ad <- Map(function(seq_data,cluster) {
      assign_somy(seq_data, cluster, uq=uq, lq=lq, somyl=somyl, somyu=somyu)
    }, peaks[, .SD, .SDcols = patterns("cell-")], clusters_pruned)
    print("Successfully assigned gain-loss")
    if(is.null(title_karyo)){
      title_karyo <- basename(outdir)
    }
    plot_karyo(somies_ad = somies_ad, outdir = outdir, peaks = peaks, uq, lq, somyl, somyu, title_karyo)
    print("Successfully plotted karyogram")
  }
}
