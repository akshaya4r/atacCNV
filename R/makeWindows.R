makeWindows <- function(genome, blacklist, windowSize, slidingSize = 2e6, exclude = c("chrM","chrX","chrY")){
  #chromSizes <- GRanges(names(seqlengths(genome)), IRanges(1, seqlengths(genome)))
  #chromSizes <- GenomeInfoDb::keepStandardChromosomes(chromSizes, pruning.mode = "coarse")
  #windows <- slidingWindows(x = chromSizes, width = windowSize, step = slidingSize) %>% unlist %>% .[which(width(.)==windowSize),]
  windows <- tileGenome(seqlengths = seqlengths(genome), tilewidth = windowSize, cut.last.tile.in.chrom = TRUE)
  windows <- GenomeInfoDb::keepStandardChromosomes(windows, pruning.mode = "coarse")
  # windows <- windows[seqnames(windows) %ni% exclude]
  windows <- dropSeqlevels(windows, exclude, pruning.mode = 'coarse')
  mcols(windows)$wSeq <- as.character(seqnames(windows))
  mcols(windows)$wStart <- BiocGenerics::start(windows)
  mcols(windows)$wEnd <- BiocGenerics::end(windows)
  message("Subtracting Blacklist...")
  # windowsBL <- lapply(seq_along(windows), function(x){
  #   if(x %% 100 == 0){
  #     message(sprintf("%s of %s", x, length(windows)))
  #   }
  #   gr <- GenomicRanges::setdiff(windows[x,], blacklist)
  #   mcols(gr) <- mcols(windows[x,])
  #   return(gr)
  # })
  # names(windowsBL) <- paste0("w",seq_along(windowsBL))
  # windowsBL <- unlist(GenomicRangesList(windowsBL), use.names = TRUE)
  # mcols(windowsBL)$name <- names(windowsBL)
  overlaps <- findOverlaps(windows, blacklist)
  idx <- setdiff(1:length(windows), S4Vectors::queryHits(overlaps))
  windowsBL <- windows[idx]
  names(windowsBL) <- paste0("w",seq_along(windowsBL))
  mcols(windowsBL)$name <- names(windowsBL)
  message("Adding Nucleotide Information...")
  windowSplit <- split(windowsBL, as.character(seqnames(windowsBL)))
  windowNuc <- lapply(seq_along(windowSplit), function(x){
    message(sprintf("%s of %s", x, length(windowSplit)))
    # chrSeq <- Biostrings::getSeq(genome,chromSizes[which(seqnames(chromSizes)==names(windowSplit)[x])])
    chrSeq <- Biostrings::getSeq(genome,names(windowSplit)[x])
    grx <- windowSplit[[x]]
    aFreq <- Biostrings::alphabetFrequency(Biostrings::Views(chrSeq, ranges(grx)))
    mcols(grx)$GC <- rowSums(aFreq[, c("G","C")]) / rowSums(aFreq)
    mcols(grx)$AT <- rowSums(aFreq[, c("A","T")]) / rowSums(aFreq)
    tn5motif1 <- Biostrings::DNAString("GSSCTGGGS")
    tn5motif2 <- Biostrings::reverseComplement(tn5motif1)
    tn5bias1 <- Biostrings::vcountPattern(tn5motif1, Biostrings::Views(chrSeq, ranges(grx)), fixed=FALSE)
    tn5bias2 <- Biostrings::vcountPattern(tn5motif2, Biostrings::Views(chrSeq, ranges(grx)), fixed=FALSE)
    mcols(grx)$tn5bias <- tn5bias1 + tn5bias2
    # grx$tn5bias[which(grx$tn5bias > quantile(grx$tn5bias, 0.90))] <- 0
    return(grx)
  }) %>% GRangesList %>% unlist %>% sortSeqlevels %>% sort
  windowNuc$N <- 1 - (windowNuc$GC + windowNuc$AT)
  windowNuc
}
