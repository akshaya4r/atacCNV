#' @export
#' @param reads bamfile or fragments file of sequencing information
#' @param windows Binned genome

generateCountMatrix <- function(reads, windows, by=NULL, minFrags=NULL){

  #Keep only regions in filtered chromosomes
  windows   <- GenomeInfoDb::keepStandardChromosomes(windows, pruning.mode = "coarse")

  #Count Insertions in windows
  message("Getting Counts...")
  counts <- countInsertions(reads, windows, by="barcode", minFrags=minFrags)

  #Keep only regions with less than 0.1% N
  keep <- which(windows$N < 0.001)
  windowSummary <- windows[keep,]
  countSummary <- counts[keep,]

  se <- SummarizedExperiment(assays = list(counts = as.matrix(countSummary)), rowRanges = windowSummary)
  colnames(se) <- colnames(counts)

  return(se)
}

# generateCountMatrix_accu <- function(reads, windows, by=NULL, minFrags=NULL){
#
#   #Keep only regions in filtered chromosomes
#   windows   <- GenomeInfoDb::keepStandardChromosomes(windows, pruning.mode = "coarse")
#
#   #Count Insertions in windows
#   message("Getting Counts...")
#   counts <- countInsertions(reads, windows, by="barcode", minFrags=minFrags)
#   message("Summarizing...")
#   windowSummary <- GenomicRangesList()
#   countSummary <- matrix(nrow=length(unique(windows$name)), ncol = ncol(counts))
#   for(x in seq_along(unique(mcols(windows)$name))){
#     if(x %% 100 == 0){
#       message(sprintf("%s of %s", x, length(unique(mcols(windows)$name))))
#     }
#     idx <- which(mcols(windows)$name == unique(mcols(windows)$name)[x])
#     wx <- windows[idx,]
#     wo <- GRanges(mcols(wx)$wSeq , ranges = IRanges(mcols(wx)$wStart, mcols(wx)$wEnd))[1,]
#     mcols(wo)$name <- mcols(wx)$name[1]
#     mcols(wo)$effectiveLength <- sum(width(wx))
#     mcols(wo)$percentEffectiveLength <- 100*sum(width(wx))/width(wo)
#     mcols(wo)$GC <- sum(mcols(wx)$GC * width(wx))/width(wo)
#     mcols(wo)$AT <- sum(mcols(wx)$AT * width(wx))/width(wo)
#     mcols(wo)$N <- sum(mcols(wx)$N * width(wx))/width(wo)
#     countSummary[x,] <- Matrix::colSums(counts[idx,,drop=FALSE])
#     windowSummary[[x]] <- wo
#   }
#   windowSummary <- unlist(windowSummary)
#
#   # Keep only regions with less than 0.1% N
#   keep <- which(windowSummary$N < 0.001)
#   windowSummary <- windowSummary[keep,]
#   countSummary <- countSummary[keep,]
#   # keep <- which(windows$N < 0.001)
#   # windowSummary <- windows[keep,]
#   # countSummary <- counts[keep,]
#   # countSummary <- countSummary[keep,]
#
#   se <- SummarizedExperiment(assays = list(counts = as.matrix(countSummary)), rowRanges = windowSummary)
#   colnames(se) <- colnames(counts)
#
#   return(se)
# }
