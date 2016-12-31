## Time-stamp: <getElbows.R zma 2010-01-17 12:05>

#' Given a decreasingly sorted vector, return the given number of elbows
#'
#' Reference:
#'   Zhu, Mu and Ghodsi, Ali (2006), "Automatic dimensionality selection from
#'   the scree plot via the use of profile likelihood", Computational
#'   Statistics & Data Analysis
#'
#' @param d: the decreasingly sorted vector (e.g. a vector of standard deviations)
#' @param n: the number of returned elbows.
#' @param threshold: either FALSE or a number. If threshold is a number, then all
#'   the elements in d that are not larger than the threshold will be ignored.
#' @return a vector of length n.
#' @author Youngser Park <youngser@@jhu.edu>
getElbows <- function(d, n = 3, threshold = FALSE)
{
  if (is.unsorted(-d))
    stop("d must be sorted decreasingly!")

  if (!is.logical(threshold))
    d <- d[d > threshold]
    
  p <- length(d)
  if (p == 0)
    stop(paste("d must have elements that are larger than the threshold ",
               threshold), "!", sep="")

  lq <- rep(0.0, p)                     # log likelihood, function of q
  for (q in 1:p) {
    mu1 <- mean(d[1:q])
    mu2 <- mean(d[-(1:q)])              # = NaN when q = p
    sigma2 <- (sum((d[1:q] - mu1)^2) + sum((d[-(1:q)] - mu2)^2)) /
      (p - 1 - (q < p))
    lq[q] <- sum( dnorm(  d[1:q ], mu1, sqrt(sigma2), log=TRUE) ) +
      sum( dnorm(d[-(1:q)], mu2, sqrt(sigma2), log=TRUE) )
  }

  q <- which.max(lq)
  if (n > 1 && q < p) {
    return(c(q, q + getElbows(d[(q+1):p], n-1)))
  } else {
    return(q)
  }
}
