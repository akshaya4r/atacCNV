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

message("--------------- Calling CNV from scATAC-seq ---------------")

atacCNV <- function(input, outdir, blacklist, windowSize){

  if(!file.exists(file.path(outdir,"count_summary.rds"))) {
    blacklist <- read_bed(blacklist)
    windows <- makeWindows(genome = BSgenome.Hsapiens.UCSC.hg38, blacklist = blacklist, windowSize)
    print(outdir)

    if(file_test("-d", input)){
      bamfiles <- Rsamtools::BamFileList(list.files(input, pattern = ".bam$", full.names = TRUE), yieldSize=100000)
      print(bamfiles)
      counts <- generateCountMatrix(windows, bamfiles, remove = c("chrM","chrX","chrY"))
    }
    else if(file_test("-f", input)){
      file_fragments <- fread(input)
      colnames(file_fragments) <- c('seqnames','start','end','barcode','pcr')
      fragments <- as_granges(file_fragments)
    }

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
