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
	# If we need to compute more right singular vectors,
	# we compute right singular vectors first.
	if (nv > nu)
		comp.right <- TRUE
	# If we need to compute more left singular vectors,
	# we compute left singular vectors first.
	else if (nv < nu)
		comp.right <- FALSE
	# If we compute the same number of left and right singular vectors,
	# we choose to compute on the form that reduces the length of vectors.
	else if (nrow(x) > ncol(x) && nv > 0)
		comp.right <- TRUE
	else
		comp.right <- FALSE

	if (comp.right) {
		n <- ncol(x)
		nev <- nv
		multiply <- function(mat, extra) tx %*% (x %*% mat)
	}
	else {
		n <- nrow(x)
		nev <- nu
		multiply <- function(mat, extra) x %*% (tx %*% mat)
	}
	res <- fm.eigen(multiply, sym=TRUE,
					options=list(n=n, nev=nev, tol=tol, which="LM"))

	# After we compute the other singular vectors, we need to rescale
	# these singular vectors because they aren't orthonormal.
	rescale <- function(x) {
		scal <- sqrt(colSums(x * x))
		x <- fm.mapply.row(x, scal, fm.bo.div)
		x <- fm.materialize(x)
	}
	if (comp.right) {
		right <- res$vecs
		left <- NULL
		if (nu > 0) {
			left <- x %*% right
			if (ncol(left) > nu)
				left <- rescale(fm.get.cols(left, 1:nu))
		}
	}
	else {
		left <- res$vecs
		right <- NULL
		if (nv > 0) {
			right <- t(x) %*% left
			if (ncol(right) > nv)
				right <- rescale(fm.get.cols(right, 1:nv))
		}
	}
	list(d=sqrt(res$vals), u=left, v=right)
}
