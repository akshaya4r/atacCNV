#' @export

assign_somy <- function(seq_data, cluster, CNgrid.start=1.5){
  counts.normal <- seq_data
  cnmean <- sapply(split(counts.normal,cluster), function(x) {
    qus <- stats::quantile(x, c(0.01, 0.99))
    y <- x[x >= qus[1] & x <= qus[2]]
    y <- y[which(y > stats::quantile(y, 0.75))]
    if(sum(y) == 0 | length(y) == 0)
      y <- x
    mean(y)
  })
  cnmean <- cnmean / min(cnmean[cnmean>0])
  # cnmean <- round(plyr::round_any(plyr::round_any(cnmean, 0.25), 0.5))
  counts.normal.mean <- cnmean[as.character(cluster)]
  CNgrid <- seq(CNgrid.start,6,by=0.01) # --> Determine copy number
  outerRaw <- counts.normal.mean %o% CNgrid
  outerDiff <- (outerRaw - round(outerRaw)) ^ 2
  sumOfSquares <- colSums(outerDiff,na.rm=FALSE,dims=1)
  CN <- CNgrid[order(sumOfSquares)][1]
  CN.states <- round(counts.normal.mean * CN, digits = 0)
  return(CN.states)
}
