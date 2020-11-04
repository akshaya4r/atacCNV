#' @export

getbp <- function(seq_data, minsize, k, test='AD'){
  bp <- vector()
  distbp <- vector()
  for(iter in 1:k) {
    seq_k_data <- splitAt(seq_data, (bp))
    dist_vect <- lapply(seq_k_data, function(x) {seq_dist_ad(x, minsize, test)})
    if(length(bp)==0){
      add_to_bp <- 0
    } else{
      add_to_bp <- c(0, sapply(seq_k_data[-length(seq_k_data)], length))
    }
    bp_per_seq <- sapply(dist_vect, function(x) { which.max(x)*minsize })
    # bp_neighbors <-  as.vector(sapply(bp, function(x){seq(x-5,x+5)}))
    # bp_per_seq[((bp_per_seq+add_to_bp) %in% bp_neighbors)] <- NA
    bp_per_seq[which(bp_per_seq > sapply(seq_k_data, length))] <- sapply(seq_k_data, length)[which(bp_per_seq > sapply(seq_k_data, length))]
    sig_per_seq <- sapply(dist_vect, max)>10
    bp_per_seq <- bp_per_seq + add_to_bp
    if(sum(sig_per_seq!=0)){
      bp <- append(bp, bp_per_seq[sig_per_seq])
      distbp <- append(distbp, sapply(dist_vect[sig_per_seq], max))
    }
  }
  res.df <- data.frame(bp, distbp)
  res.df <- res.df[order(-res.df$distbp),]
  res <- list()
  res[['bp']] <- res.df$bp
  res[['ad.distance']] <- res.df$distbp
  # res[['estimates']] <- unique(res$bp[1:(elbow_finder(res$ad.distance)-1)])
  bp_to_cluster <- sort(c(1,length(seq_data)+1,res$bp))
  res[['cluster']] <- rep(1:length(diff(bp_to_cluster)),diff(bp_to_cluster))
  # res[['somy']] <- assign_somy(seq_data,res$cluster)
  return(res)
}
