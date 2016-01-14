rank_matrices <- function(graphs, normalize=FALSE) {
  d <- dim(graphs)
  rg <- array(rep(NaN, d[1]*d[2]*d[3]), d)
  for (i in 1:d[3]) {
    rg[,,i] <- array(rank(graphs[,,i], ties.method="average"), c(d[1], d[2]))
    if (normalize) {
      rg[,,i] <- ( rg[,,i] - min(rg[,,i]) ) / (max(rg[,,i]) - min(rg[,,i]))
    }
  }
  return(rg)
}
