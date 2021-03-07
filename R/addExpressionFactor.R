#' @import GenomicFeatures
#' @export

addExpressionFactor <- function(bins, gene.annotation=NULL) {
  UseMethod("addExpressionFactor")
}

# Looking at genes
addExpressionFactor.GRanges <- function(bins, gene.annotation=NULL) {
  txdb <- getFromNamespace(gene.annotation, ns=gene.annotation)
  seqlevelsStyle(txdb) <- seqlevelsStyle(bins)[1]
  genes <- sort(keepStandardChromosomes(genes(txdb), pruning.mode = 'coarse'))
  bins$genecoverage <- countOverlaps(bins,genes)
  # bins[is.na(bins$genecoverage)]$genecoverage <- 0
  bins
}
# addExpressionFactor.GRanges <- function(bins, gene.annotation=NULL) {
#   txdb <- getFromNamespace(gene.annotation, ns=gene.annotation)
#   seqlevelsStyle(txdb) <- seqlevelsStyle(bins)[1]
#   genes <- sort(keepStandardChromosomes(genes(txdb), pruning.mode = 'coarse'))
#   genes <- reduce(genes, ignore.strand=TRUE)
#   hits <- findOverlaps(bins, genes)
#   overlaps <- as.data.table(pintersect(bins[queryHits(hits)], genes[subjectHits(hits)]))
#   overlaps <- overlaps[ , .(genecoverage=sum(width)), by=name]
#   print(overlaps)
#   bins <- sort(as_granges(merge(bins, overlaps, all=TRUE)))
#   bins[is.na(bins$genecoverage)]$genecoverage <- 0
#   bins
# }

addExpressionFactor.list <- function(bins, gene.annotation=NULL) {
  lapply(bins, addExpressionFactor, txdb)
}

addExpressionFactor.default <- function(bins, gene.annotation=NULL) {
  stop("Do not know how to get expression factor for type ", sQuote(class(bins)))
}
