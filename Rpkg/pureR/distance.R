compute_distance <- function(graphs, normx='F') {
  S <- dim(graphs)[3]
  dist <- matrix(rep(0, S*S), ncol=S)
  for (i in 1:dim(graphs)[3]) {
    for (j in i:dim(graphs)[3]) {
      dist[i,j] <- norm(graphs[,,i]-graphs[,,j], normx)
    }
  }
  dist <- dist + t(dist)
  return(dist)
}
