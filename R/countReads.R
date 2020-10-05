countInsertions <- function(windows, bamfiles, by = "RG"){
  # seqlevelsStyle(windows) <- 'NCBI'
  counts_bam <- summarizeOverlaps(windows, bamfiles, singleEnd=TRUE, fragments=FALSE, mode='Union', param = ScanBamParam(mapqFilter=10), ignore.strand = TRUE)
  sparseM <- Matrix(assays(counts_bam)$counts, sparse = TRUE)
  frip <- 1
  total <- colSums(sparseM)
  out <- list(counts = sparseM, frip = frip, total = total)
  return(out)
}
