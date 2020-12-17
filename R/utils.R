splitAt <- function(x, pos) unname(split(x, cumsum(seq_along(x) %in% pos)))

# "%ni%" <- function(){ Negate("%in%") }

qc.spikiness <- function(counts) {
  if (is.null(counts)) {
    return(NA)
  }
  counts <- as.vector(counts)
  sum.counts <- sum(counts)
  spikiness <- sum(abs(diff(counts))) / sum.counts
  return(spikiness)
}

qc.entropy <- function(counts) {
  if (is.null(counts)) {
    return(NA)
  }
  counts <- as.vector(counts)
  total.counts <- sum(counts)
  n <- counts/total.counts
  entropy <- -sum( n * log(n) , na.rm=TRUE)
  return(entropy)
}

qc.sos <- function(counts, somies) {
  sum(counts - somies) ^ 2
}
