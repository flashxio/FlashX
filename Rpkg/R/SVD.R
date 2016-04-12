#' Compute the singular-value decomposition of a large sparse matrix.
#'
#' @param x a FlashMatrixR object
#' @param nu the number of left singluar vectors to be computed.
#' @param nv the number of right singluar vectors to be computed.
#' @param tol Stopping criterion: the relative accuracy of the Ritz value
#' is considered acceptable if its error is less than 'tol' times its estimated
#' value. If this is set to zero then machine precision is used.
#' @return Returns a list with three entries
#' \itemize{
#'   \item{d}{ max(nu, nv) approximate singular values}
#'   \item{u}{ nu approximate left singular vectors (only when right_only=FALSE)}
#'   \item{v}{ nv approximate right singular vectors}
#' }
#' @author Da Zheng <dzheng5@@jhu.edu>
fm.svd <- function(x, nu, nv, tol=1e-8)
{
	stopifnot(class(x) == "fm")
	tx <- t(x)
	comp.right <- FALSE
	if (nrow(x) > ncol(x))
		comp.right <- TRUE

	nev <- max(nu, nv)
	if (comp.right) {
		n <- ncol(x)
		x.prod <- fm.conv.FM2R(tx %*% x)
		multiply <- function(vec, extra) x.prod %*% vec
	}
	else {
		n <- nrow(x)
		x.prod <- fm.conv.FM2R(x %*% tx)
		multiply <- function(vec, extra) x.prod %*% vec
	}
	# If it's a very small matrix, we can compute its eigenvalues directly.
	# Or if we need to compute many eigenvalues, we probably should also
	# compute its eigenvalues directly.
	if (nrow(x.prod) < 100 || nev >= nrow(x.prod) / 2) {
		res <- eigen(x.prod, TRUE, FALSE)
		res$values <- res$values[nev]
		nev <- nrow(x.prod)
	}
	else
		res <- arpack(multiply, sym=TRUE,
					  options=list(n=n, nev=nev, tol=tol, ncv=max(nev * 2, 5), which="LM"))

	# After we compute the other singular vectors, we need to rescale
	# these singular vectors because they aren't orthonormal.
	rescale <- function(x) {
		fm.set.materialize.level(x, 2)
		scal <- sqrt(colSums(x * x))
		x <- fm.mapply.row(x, scal, fm.bo.div)
#		x <- fm.materialize(x)
	}
	if (comp.right) {
		right <- NULL
		if (nv > 0) {
			if (nv < nev)
				right <- fm.conv.R2FM(res$vectors[,1:nv])
			else
				right <- fm.conv.R2FM(res$vectors)
		}
		left <- NULL
		if (nu > 0) {
			left <- x %*% res$vectors[,1:nu]
			left <- rescale(left)
		}
	}
	else {
		left <- NULL
		if (nu > 0) {
			if (nu < nev)
				left <- fm.conv.R2FM(res$vectors[,1:nu])
			else
				left <- fm.conv.R2FM(res$vectors)
		}
		right <- NULL
		if (nv > 0) {
			right <- t(x) %*% res$vectors[,1:nv]
			right <- rescale(right)
		}
	}
	list(d=sqrt(res$values), u=left, v=right)
}
