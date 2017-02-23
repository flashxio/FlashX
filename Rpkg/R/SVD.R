# Copyright 2016 Open Connectome Project (http://openconnecto.me)
# Written by Da Zheng (zhengda1936@gmail.com)
#
# This file is part of FlashR.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

#' Compute the singular-value decomposition on a large matrix.
#'
#' The difference between \code{svd} and \code{fm.svd} is that \code{fm.svd}
#' allows a user-specified tol. \code{svd} computes eigenvalues in machine
#' precision.
#'
#' @param x a FlashR matrix
#' @param nu the number of left singluar vectors to be computed.
#' @param nv the number of right singluar vectors to be computed.
#' @param tol Stopping criterion: the relative accuracy of the Ritz value
#' is considered acceptable if its error is less than 'tol' times its estimated
#' value. If this is set to zero then machine precision is used.
#' @return Returns a list with three entries
#'   \item{d}{ max(nu, nv) approximate singular values}
#'   \item{u}{ nu approximate left singular vectors (only when right_only=FALSE)}
#'   \item{v}{ nv approximate right singular vectors}
#' @author Da Zheng <dzheng5@@jhu.edu>
#' @name svd
#'
#' @examples
#' mat <- fm.runif.matrix(1000, 100)
#' res <- fm.svd(mat, nu=10, nv=0)
#' res <- svd(mat, nu=10, nv=0)
NULL

#' @rdname svd
fm.svd <- function(x, nu=min(n, p), nv=min(n, p), tol=1e-8)
{
	stopifnot(class(x) == "fm")
	x <- as.double(x)
	n <- dim(x)[1]
	p <- dim(x)[2]
	tx <- t(x)
	comp.right <- FALSE
	if (nrow(x) > ncol(x))
		comp.right <- TRUE

	nev <- max(nu, nv)
	x.prod <- NULL
	if (fm.is.sparse(x) && comp.right) {
		size <- ncol(x)
		multiply <- function(vec, extra) t(x) %*% (x %*% vec)
	}
	else if (fm.is.sparse(x)) {
		size <- nrow(x)
		multiply <- function(vec, extra) x %*% (t(x) %*% vec)
	}
	else if (comp.right) {
		size <- ncol(x)
		x.prod <- tx %*% x
		multiply <- function(vec, extra) x.prod %*% vec
	}
	else {
		size <- nrow(x)
		x.prod <- x %*% tx
		multiply <- function(vec, extra) x.prod %*% vec
	}
	# If it's a very small matrix, we can compute its eigenvalues directly.
	# Or if we need to compute many eigenvalues, we probably should also
	# compute its eigenvalues directly.
	if (!is.null(x.prod) && (size < 100 || nev >= size / 2)) {
		if (!is.matrix(x.prod))
			x.prod <- as.matrix(x.prod)
		res <- eigen(x.prod, TRUE, FALSE)
		res$values <- res$values[1:nev]
		nev <- nrow(x.prod)
		res$vectors <- fm.as.matrix(res$vectors)
	}
	else
		res <- fm.eigen(multiply, nev, size, which="LM", sym=TRUE,
					  options=list(tol=tol, ncv=max(nev * 2, 5)))
	if (fm.is.vector(res$vectors))
		res$vectors <- fm.as.matrix(res$vectors)

	# After we compute the other singular vectors, we need to rescale
	# these singular vectors because they aren't orthonormal.
	rescale <- function(x) {
		if (fm.is.vector(x))
			x <- fm.as.matrix(x)
		fm.set.cached(x, TRUE)
		scal <- sqrt(colSums(x * x))
		x <- fm.mapply.row(x, scal, fm.bo.div)
	}
	if (comp.right) {
		right <- NULL
		if (nv > 0) {
			if (nv < nev)
				right <- res$vectors[,1:nv]
			else
				right <- res$vectors
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
				left <- res$vectors[,1:nu]
			else
				left <- res$vectors
		}
		right <- NULL
		if (nv > 0) {
			right <- t(x) %*% res$vectors[,1:nv]
			right <- rescale(right)
		}
	}
	# If an eigenvalue is very small (close to the machine precision), it's
	# possible that the eigenvalue is negative but very close to 0.
	vals <- ifelse(res$values > 0, res$values, 0)
	list(d=sqrt(vals), u=left, v=right, options=res$options)
}

#' @rdname svd
setMethod("svd", signature(x = "fm"), function(x, nu=min(n, p), nv=min(n, p), LINPACK) {
		  x <- fm.as.matrix(x)
		  if (any(!is.finite(x)))
			  stop("infinite or missing values in 'x'")
		  dx <- dim(x)
		  n <- dx[1L]
		  p <- dx[2L]
		  if (!n || !p)
			  stop("a dimension is zero")
		  fm.res <- fm.svd(x, nu, nv, tol=.Machine$double.eps)
		  res <- list(d = fm.res$d)
		  if (nu)
			  res$u <- fm.res$u
		  if (nv)
			  res$v <- fm.res$v
		  res
})

setMethod("prcomp", signature(x = "fm"), function(x, retx=TRUE, center=TRUE,
												  scale.=FALSE, tol=NULL) {
	scale.x <- scale(x, center, scale.)
	res <- fm.svd(scale.x, nu=0, tol=.Machine$double.eps)
	if (!is.null(tol)) {
		idxs <- which(res$d > tol)
		rotation <- res$v[, idxs]
	}
	else
		rotation <- res$v
	if (retx)
		x <- scale.x %*% rotation
	else
		x <- NULL
	list(sdev=res$d / sqrt(nrow(scale.x) - 1), rotation=rotation, x=x)
})
