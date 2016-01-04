procrustes <- function(X,Y, type = "I", normFlag=FALSE) {
  if(type == "C"){
    X <- X/norm(X, type = "F")*sqrt(nrow(X))
    Y <- Y/norm(Y, type = "F")*sqrt(nrow(Y))
  }
  if(type == "D"){
    tX <- rowSums(X^2)
    tX[tX <= 1e-15] <- 1
    tY <- rowSums(Y^2)
    tY[tY <= 1e-15] <- 1
    X <- X/sqrt(tX)
    Y <- Y/sqrt(tY)
  }

  tmp <- t(X) %*% Y
  tmp.svd <- svd(tmp)
  W <- tmp.svd$u %*% t(tmp.svd$v)
  err <- norm(X%*%W - Y, type = "F")
  if (normFlag) err <- err / (computeCX(X) + computeCX(Y))

  return(list(error = err, W = W))
}
