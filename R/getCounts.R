#' @export
#' @param reads bamfile or fragments file
#' @param windows binned genome as generated by the makeWindows function
#' @param by barcode information
#' @param minFrags Minimum number of fragments a barcode must contain to be counted as a cell
#'
countInsertions <- function(reads, windows, by=NULL, minFrags=NULL) {
  UseMethod("countInsertions")
}

countInsertions.BamFileList <- function(reads, windows, by=NULL, minFrags=NULL){
  # seqlevelsStyle(windows) <- 'NCBI'
  message("Counting reads from bam files .. ")
  counts_bam <- GenomicAlignments::summarizeOverlaps(windows, reads, singleEnd=TRUE, fragments=FALSE, mode='Union', param = Rsamtools::ScanBamParam(mapqFilter=10), ignore.strand = TRUE)
  sparseM <- Matrix(assays(counts_bam)$counts, sparse = TRUE)
  # frip <- 1
  # total <- colSums(sparseM)
  # out <- list(counts = sparseM, frip = frip, total = total)
  return(sparseM)
}

countInsertions.GRanges <- function(reads, windows, by = "barcode", minFrags = 20000){
  message("Counting reads from fragments/bed file .. ")
  tabRG <- table(mcols(reads)[[by]])
  keep <- names(tabRG)[which(tabRG >= minFrags)]
  reads <- reads[mcols(reads)[[by]] %in% keep,]
  reads <- sort(sortSeqlevels(reads))
  overlapDF <- DataFrame(findOverlaps(windows, reads, ignore.strand = TRUE, maxgap=-1L, minoverlap=0L, type = "any"))
  overlapDF$name <- mcols(reads)[overlapDF[, 2], by]
  overlapTDF <- transform(overlapDF, id = match(name, unique(name)))
  sparseM <- Matrix::sparseMatrix(
    i = overlapTDF[, 1],
    j = overlapTDF[, 4],
    x = rep(1, nrow(overlapTDF)),
    dims = c(length(windows), length(unique(overlapDF$name))))
  colnames(sparseM) <- unique(overlapDF$name)
  return(sparseM)
}
