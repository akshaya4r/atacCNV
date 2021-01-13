#' @export

threshold_dist_values <- function(result.dt) {
  result.dt$zscores <- scale(result.dt$per_chr_dist, center = TRUE, scale = TRUE)
  result.dt <- result.dt[zscores>0,]
  result.dt$zscores <- NULL
  return(result.dt)
}

# assign_somy <- function(seq_data, cluster, CNgrid.start=1.5){
#   counts.normal <- seq_data
#   cnmean <- sapply(split(counts.normal,cluster), function(x) {
#     qus <- stats::quantile(x, c(0.01, 0.99))
#     y <- x[x >= qus[1] & x <= qus[2]]
#     y <- y[which(y > stats::quantile(y, 0.75))]
#     if(sum(y) == 0 | length(y) == 0)
#       y <- x
#     mean(y)
#   })
#   cnmean <- cnmean / min(cnmean[cnmean>0])
#   # cnmean <- round(plyr::round_any(plyr::round_any(cnmean, 0.25), 0.5))
#   counts.normal.mean <- cnmean[as.character(cluster)]
#   CNgrid <- seq(CNgrid.start,6,by=0.01) # --> Determine copy number
#   outerRaw <- counts.normal.mean %o% CNgrid
#   outerDiff <- (outerRaw - round(outerRaw)) ^ 2
#   sumOfSquares <- colSums(outerDiff,na.rm=FALSE,dims=1)
#   CN <- CNgrid[order(sumOfSquares)][1]
#   CN.states <- round(counts.normal.mean * CN, digits = 0)
#   return(CN.states)
# }

assign_somy <- function(seq_data, cluster, CNgrid.start=1.5, uq=0.8, lq=0.5, somyl=0.2, somyu=0.8){
  counts.normal <- seq_data
  qus_global <- quantile(seq_data, c(0.01, 0.98))
  # counts.normal <- seq_data[seq_data >= qus[1] & seq_data <= qus[2]]
  # cluster_quantile <- cluster[seq_data >= qus[1] & seq_data <= qus[2]]
  cnmean <- sapply(split(counts.normal,cluster), function(x) {
    qus <- quantile(x, c(lq, uq))
    # qus <- quantile(x, c(0.50, 0.99))
    y <- x[x >= qus[1] & x <= qus[2] & x >= qus_global[1] & x <= qus_global[2]]
    # y <- y[which(y > quantile(y, 0.70))]
    if(sum(y) == 0 | length(y) == 0)
      y <- x
    mean(y)
  })
  raw_cnmean <- cnmean / mean(cnmean)
  # raw_cnmean[raw_cnmean>quantile(raw_cnmean, 0.80) | raw_cnmean<quantile(raw_cnmean, 0.20)] <- NA
  # cnmean[cnmean>quantile(cnmean, 0.80)] <- quantile(cnmean, 0.80)
  # cnmean[cnmean>quantile(cnmean, 0.75) | cnmean<quantile(cnmean, 0.25)] <- mean(cnmean)
  cnmean[cnmean>quantile(cnmean, somyu) | cnmean<quantile(cnmean, somyl)] <- NA
  # cnmean <- cnmean / min(cnmean[cnmean>0])
  # cnmean <- runmed(cnmean, k=3)
  # cooksd <- cooks.distance(lm(cnmean~1))
  # cnmean[(cooksd > (1/length(cnmean)))] <- mean(cnmean)
  cnmean <- cnmean / mean(cnmean, na.rm=TRUE)
  # cnmean <- round(plyr::round_any(plyr::round_any(cnmean, 0.25), 0.5))
  counts.normal.mean <- cnmean[as.character(cluster)]
  CNgrid <- seq(CNgrid.start,6,by=0.01) # --> Determine copy number
  outerRaw <- counts.normal.mean %o% CNgrid
  outerDiff <- (outerRaw - round(outerRaw)) ^ 2
  sumOfSquares <- colSums(outerDiff,na.rm=TRUE,dims=1)
  CN <- CNgrid[order(sumOfSquares)][1]
  # CN
  CN.states <- round(raw_cnmean[as.character(cluster)] * CN)#, digits = 0)
  cnmean_intsomy <- round(raw_cnmean*CN)
  # cnmean_intsomy
  # return(list(raw_cnmean, cnmean, CN, cnmean_intsomy))
  return(CN.states)
}

assign_gainloss <- function(seq_data, cluster, CNgrid.start=1.5, uq=0.8, lq=0.1) {
  counts.normal <- seq_data
  qus_global <- quantile(seq_data, c(0.01, 0.98))
  # counts.normal <- seq_data[seq_data >= qus[1] & seq_data <= qus[2]]
  # cluster_quantile <- cluster[seq_data >= qus[1] & seq_data <= qus[2]]
  cnmean <- sapply(split(counts.normal,cluster), function(x) {
    qus <- quantile(x, c(lq, uq))
    # qus <- quantile(x, c(0.50, 0.99))
    y <- x[x >= qus[1] & x <= qus[2] & x >= qus_global[1] & x <= qus_global[2]]
    # y <- y[which(y > quantile(y, 0.70))]
    if(sum(y) == 0 | length(y) == 0)
      y <- x
    mean(y)
  })
  raw_cnmean <- cnmean / mean(cnmean)
  CN.states <- round(raw_cnmean[as.character(cluster)])
  # CN.states <- raw_cnmean[as.character(cluster)]
  return(CN.states)
}

