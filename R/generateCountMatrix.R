#' @export

generateCountMatrix <- function(reads, windows, ...){

  arguments <- list(...)
  #Keep only regions in filtered chromosomes
  windows   <- GenomeInfoDb::keepStandardChromosomes(windows, pruning.mode = "coarse")

  #Count Insertions in windows
  message("Getting Counts...")
  counts <- countInsertions(reads, windows, arguments)

  #Keep only regions with less than 0.1% N
  keep <- which(windows$N < 0.001)
  windowSummary <- windows[keep,]
  countSummary <- counts[keep,]

  se <- SummarizedExperiment(assays = list(counts = as.matrix(countSummary)), rowRanges = windowSummary)
  colnames(se) <- colnames(counts)

  return(se)
}

# generateCountMatrix <- function(windows, reads, remove = c("chrM","chrX","chrY")){
#
#   #Keep only regions in filtered chromosomes
#   windows   <- GenomeInfoDb::keepStandardChromosomes(windows, pruning.mode = "coarse")
#   # fragments <- GenomeInfoDb::keepStandardChromosomes(fragments, pruning.mode = "coarse")
#   # windows <- windows[seqnames(windows) %ni% remove]
#   # fragments <- fragments[seqnames(fragments) %ni% remove]
#
#   #Count Insertions in windows
#   message("Getting Counts...")
#   counts <- countInsertions(windows, reads, by = "RG")[[1]]
#   message("Summarizing...")
#   # windowSummary <- GenomicRangesList()
#   # countSummary <- matrix(nrow=length(unique(windows$name)), ncol = ncol(counts))
#   # for(x in seq_along(unique(mcols(windows)$name))){
#   #   if(x %% 100 == 0){
#   #     message(sprintf("%s of %s", x, length(unique(mcols(windows)$name))))
#   #   }
#   #   idx <- which(mcols(windows)$name == unique(mcols(windows)$name)[x])
#   #   wx <- windows[idx,]
#   #   wo <- GRanges(mcols(wx)$wSeq , ranges = IRanges(mcols(wx)$wStart, mcols(wx)$wEnd))[1,]
#   #   mcols(wo)$name <- mcols(wx)$name[1]
#   #   mcols(wo)$effectiveLength <- sum(width(wx))
#   #   mcols(wo)$percentEffectiveLength <- 100*sum(width(wx))/width(wo)
#   #   mcols(wo)$GC <- sum(mcols(wx)$GC * width(wx))/width(wo)
#   #   mcols(wo)$AT <- sum(mcols(wx)$AT * width(wx))/width(wo)
#   #   mcols(wo)$N <- sum(mcols(wx)$N * width(wx))/width(wo)
#   #   mcols(wo)$tn5bias <- sum(mcols(wx)$tn5bias * width(wx))/width(wo)
#   #   countSummary[x,] <- Matrix::colSums(counts[idx,,drop=FALSE])
#   #   windowSummary[[x]] <- wo
#   # }
#   # windowSummary <- unlist(windowSummary)
#
#   #Keep only regions with less than 0.1% N
#   # keep <- which(windowSummary$N < 0.001)
#   # windowSummary <- windowSummary[keep,]
#   # countSummary <- countSummary[keep,]
#   keep <- which(windows$N < 0.001)
#   windowSummary <- windows[keep,]
#   countSummary <- counts[keep,]
#   # countSummary <- countSummary[keep,]
#
#   se <- SummarizedExperiment(assays = list(counts = as.matrix(countSummary)), rowRanges = windowSummary)
#   colnames(se) <- colnames(counts)
#
#   return(se)
# }
