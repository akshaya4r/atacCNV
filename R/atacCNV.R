# For scCAT K562, it has been aligned to the NCBI reference. Change seq level style in the countInsertion fn
<<<<<<< HEAD

# library(matrixStats)
# library(readr)
# library(edgeR)
# library(biomaRt)
# library(ggplot2)
# library(grid)
# library(cowplot)
# library(kSamples)

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
#' @import SummarizedExperiment
#' @import magrittr
#' @importFrom GenomicAlignments summarizeOverlaps

atacCNV <- function(indir, outdir, blacklist, windowSize){

  bamfiles <- Rsamtools::BamFileList(list.files(indir, pattern = ".bam$", full.names = TRUE), yieldSize=100000)
  blacklist <- read_bed(blacklist)
  print(bamfiles)
  print(outdir)

  if(!file.exists(file.path(outdir,"count_summary.rds"))) {
    windows <- makeWindows(genome = BSgenome.Hsapiens.UCSC.hg38, blacklist = blacklist, windowSize)
    counts <- generateCountMatrix(windows, bamfiles, remove = c("chrM","chrX","chrY"))
    saveRDS(counts, file.path(outdir,"count_summary.rds"))
  }

  counts <- readRDS(file.path(outdir,"count_summary.rds"))
  peaks <- as.data.table(assays(counts)$counts)
  rowinfo <- as.data.table(rowRanges(counts))
  peaks <- cbind(rowinfo, peaks)

  if(!file.exists(file.path(outdir,"counts_gc_corrected.rds"))) {
    corrected_counts <- peaks[, lapply(.SD, function(x) {
      fit <- stats::loess(x ~ peaks$GC)
      correction <- median(x) / fit$fitted
      as.integer(round(x * correction))
    }), .SDcols = patterns("bam")]
    saveRDS(corrected_counts, file.path(outdir,"counts_gc_corrected.rds"))
  }

  corrected_counts <- readRDS(file.path(outdir,"counts_gc_corrected.rds"))
  peaks <- cbind(rowinfo, corrected_counts)

  if(!file.exists(file.path(outdir,"results_gc_corrected.rds"))) {
    clusters_ad <- peaks[, lapply(.SD, function(x) {
      k <- 3 # k produces 2^k - 1 breakpoints
      minsize <- 5
      peaksperchrom <- split(x, peaks$seqnames)
      print("Calculating distance AD")
      results <- lapply(peaksperchrom, function(x2) {
        # x2 <- runmed(x2, k=11, endrule = 'median')
        getbp(x2, k = k, minsize = minsize, test='AD')
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
    }), .SDcols = patterns("bam")]
    saveRDS(clusters_ad, file.path(outdir, "results_gc_corrected.rds"))
  }
  print("Finished successfully")

  # somies_ad <- Map(function(seq_data,cluster) {
  #   assign_somy(seq_data, cluster)
  # }, peaks[, .SD, .SDcols = patterns("bam")], clusters_ad)
}
=======
library(Matrix)
library(SummarizedExperiment)
library(matrixStats)
library(readr)
library(GenomicRanges)
library(magrittr)
library(edgeR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Rsamtools)
library(GenomicAlignments)
library(biomaRt)
library(plyranges)
library(data.table)
library(ggplot2)
library(grid)
library(cowplot)
library(kSamples)

bamfiles <- BamFileList(list.files("/Users/akshaya/Work/atacCNV/data/ATAC/HCT116", pattern = ".bam$", full.names = TRUE), yieldSize=100000)
outdir <- args[2]
blacklist <- import.bed("/Users/akshaya/Work/hg38.blacklist.bed")
windowSize <- 1e6

print(args[1])
print(bamfiles)
print(outdir)

if(!file.exists(paste0(outdir,"count_summary.rds"))) {
  windows <- makeWindows(genome = BSgenome.Hsapiens.UCSC.hg38, blacklist = blacklist, windowSize)
  counts <- generateCountMatrix(windows, bamfiles, remove = c("chrM","chrX","chrY"))
  saveRDS(counts, paste0(outdir,"count_summary.rds"))
}

counts <- readRDS(paste0(outdir,"count_summary.rds"))
peaks <- as.data.table(assays(counts)$counts)
rowinfo <- as.data.table(rowRanges(counts))
peaks <- cbind(rowinfo, peaks)

if(!file.exists(paste0(outdir,"counts_gc_corrected.rds"))) {
  corrected_counts <- peaks[, lapply(.SD, function(x) {
    fit <- stats::loess(x ~ peaks$GC)
    correction <- median(x) / fit$fitted
    as.integer(round(x * correction))
  }), .SDcols = patterns("bam")]
  saveRDS(corrected_counts, paste0(outdir,"counts_gc_corrected.rds"))
}

corrected_counts <- readRDS(paste0(outdir,"counts_gc_corrected.rds"))
peaks <- cbind(rowinfo, corrected_counts)

if(!file.exists(paste0(outdir,"results_gc_corrected.rds"))) {
  clusters_ad <- peaks[, lapply(.SD, function(x) {
    k <- 3 # k produces 2^k - 1 breakpoints
    minsize <- 5
    peaksperchrom <- split(x, peaks$seqnames)
    print("Calculating distance AD")
    results <- lapply(peaksperchrom, function(x2) {
      # x2 <- runmed(x2, k=11, endrule = 'median')
      getbp(x2, k = k, minsize = minsize, test='AD')
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
  }), .SDcols = patterns("bam")]
  saveRDS(clusters_ad, paste0(outdir, "results_gc_corrected.rds"))
}
print("Finished successfully")

# somies_ad <- Map(function(seq_data,cluster) {
#   assign_somy(seq_data, cluster)
# }, peaks[, .SD, .SDcols = patterns("bam")], clusters_ad)
>>>>>>> 42fe52324aee225d722ea841f73233548074c47b
