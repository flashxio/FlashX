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

.sd.int <- function(x, na.rm)
{
	n <- length(x)
	x2 <- x * x
	test.na <- TRUE
	num.na <- 0
	if (na.rm) {
		zero <- .get.zero(typeof(x))
		in.is.na <- is.na(x)
		x <- ifelse(in.is.na, zero, x)
		x2 <- ifelse(in.is.na, zero, x2)
		test.na <- FALSE
	}
	sum.x <- fm.agg.lazy(x, fm.bo.add)
	sum.x2 <- fm.agg.lazy(x2, fm.bo.add)

	if (test.na) {
		x.is.na <- fm.agg.lazy(.is.na.only(x), fm.bo.or)
		res <- fm.materialize(sum.x, sum.x2, x.is.na)
		if (.fmV2scalar(res[[3]]))
			return(.get.na(typeof(x)))
		sums <- res[1:2]
	}
	else {
		# If we remove NA, we should calculate the number of
		# NAs in the vector.
		sum.na <- fm.agg.lazy(in.is.na, fm.bo.add)
		res <- fm.materialize(sum.x, sum.x2, sum.na)
		n <- n - .fmV2scalar(res[[3]])
		sums <- res[1:2]
	}
	sum.x <- .fmV2scalar(sums[[1]])
	sum.x2 <- .fmV2scalar(sums[[2]])
	avg <- sum.x / n
	sqrt((sum.x2 - n * avg * avg) / (n - 1))
}

#' Standard Deviation
#'
#' This function computes the standard deviation of the values in \code{x}.
#' If \code{na.rm} is \code{TRUE} then missing values are removed before
#' computation proceeds.
#'
#' @param x a FlashMatrix vector or matrix.
#' @param na.rm logical. Should missing values be removed?
#' @name sd
NULL

#' @rdname sd
setMethod("sd", "fm", .sd.int)
#' @rdname sd
setMethod("sd", "fmV", .sd.int)

.cov.int <- function(x, y=NULL, use="everything",
				   method=c("pearson", "kendall", "spearman"))
{
	orig.test.na <- .env.int$fm.test.na
	.set.test.na(FALSE)
	if (!is.null(y))
		stopifnot(nrow(x) == nrow(y))
	n <- nrow(x)
	if (is.null(y)) {
		# we need to transpose x and compute rowSum instead of computing
		# colSum on the original matrix. The reason is that fm.materialize
		# only works on a set of matrices the same the long dimension and
		# dimension size. TODO I need to change that.
		x.sum <- fm.rowSums(t(x), TRUE)
		x.prod <- fm.multiply(t(x), x, TRUE)
		ret <- fm.materialize(x.sum, x.prod)
		x.mu <- ret[[1]] / n
		x.mu <- fm.conv.FM2R(x.mu)
		ret <- (fm.conv.FM2R(ret[[2]]) - n * x.mu %*% t(x.mu)) / (n - 1)
	}
	else {
		x.sum <- fm.rowSums(t(x), TRUE)
		y.sum <- fm.rowSums(t(y), TRUE)
		xy.prod <- fm.multiply(t(x), y, TRUE)
		ret <- fm.materialize(x.sum, y.sum, xy.prod)
		x.mu <- ret[[1]] / n
		y.mu <- ret[[2]] / n
		x.mu <- fm.conv.FM2R(x.mu)
		y.mu <- fm.conv.FM2R(y.mu)
		ret <- (fm.conv.FM2R(ret[[3]]) - n * x.mu %*% t(y.mu)) / (n - 1)
	}
	.set.test.na(orig.test.na)
	ret
}

.cor.int <- function(x, y=NULL, use="everything",
				   method=c("pearson", "kendall", "spearman"))
{
	orig.test.na <- .env.int$fm.test.na
	.set.test.na(FALSE)
	if (!is.null(y))
		stopifnot(nrow(x) == nrow(y))
	n <- nrow(x)
	if (is.null(y)) {
		x.sum <- fm.rowSums(t(x), TRUE)
		x2.sum <- fm.rowSums(t(x * x), TRUE)
		x.prod <- fm.multiply(t(x), x, TRUE)
		ret <- fm.materialize(x.sum, x2.sum, x.prod)
		x.mu <- ret[[1]] / n
		x.sd <- sqrt((ret[[2]] - n * x.mu * x.mu) / (n - 1))
		x.mu <- fm.conv.FM2R(x.mu)
		x.sd <- fm.conv.FM2R(x.sd)
		ret <- (fm.conv.FM2R(ret[[3]]) - n * x.mu %*% t(x.mu)) / (n - 1) / (x.sd %*% t(x.sd))
	}
	else {
		x.sum <- fm.rowSums(t(x), TRUE)
		y.sum <- fm.rowSums(t(y), TRUE)
		x2.sum <- fm.rowSums(t(x * x), TRUE)
		y2.sum <- fm.rowSums(t(y * y), TRUE)
		xy.prod <- fm.multiply(t(x), y, TRUE)
		ret <- fm.materialize(x.sum, y.sum, x2.sum, y2.sum, xy.prod)
		x.mu <- ret[[1]] / n
		x.sd <- sqrt((ret[[3]] - n * x.mu * x.mu) / (n - 1))
		x.mu <- fm.conv.FM2R(x.mu)
		x.sd <- fm.conv.FM2R(x.sd)
		y.mu <- ret[[2]] / n
		y.sd <- sqrt((ret[[4]] - n * y.mu * y.mu) / (n - 1))
		y.mu <- fm.conv.FM2R(y.mu)
		y.sd <- fm.conv.FM2R(y.sd)
		ret <- (fm.conv.FM2R(ret[[5]]) - n * x.mu %*% t(y.mu)) / (n - 1) / (x.sd %*% t(y.sd))
	}
	.set.test.na(orig.test.na)
	ret
}

.cov.wt.int <- function (x, wt = rep(1/nrow(x), nrow(x)), cor = FALSE, center = TRUE,
					       method = c("unbiased", "ML"))
{
	orig.test.na <- .env.int$fm.test.na
	.set.test.na(FALSE)
	if (is.data.frame(x))
		x <- as.matrix(x)
	else if (!is.matrix(x))
		stop("'x' must be a matrix or a data frame")
	n <- nrow(x)
	if (with.wt <- !missing(wt)) {
		if (length(wt) != n)
			stop("length of 'wt' must equal the number of rows in 'x'")
		any.neg <- fm.any(wt < 0, TRUE)
		s <- fm.sum(wt, TRUE)
		ret <- fm.materialize(any.neg, s)
		any.neg <- fm.conv.FM2R(ret[[1]])
		s <- fm.conv.FM2R(ret[[2]])
		if (any.neg || s == 0)
			stop("weights must be non-negative and not all zero")
		wt <- wt/s
	}

	all.finite <- fm.all(is.finite(t(x)), TRUE)
	if (is.logical(center)) {
		center <- if (center)
#			fm.colSums(wt * x, TRUE)
			fm.rowSums(t(wt * x), TRUE)
		else 0
	}
	else {
		if (length(center) != ncol(x))
			stop("length of 'center' must equal the number of columns in 'x'")
	}
	wx <- wt * x
#	wx.cs <- fm.colSums(wx, TRUE)
	wx.cs <- fm.rowSums(t(wx), TRUE)
	x.cp <- fm.crossprod(wx, x, lazy=TRUE)
	if (fm.is.sink(center)) {
		ret <- fm.materialize(center, all.finite, wx.cs, x.cp)
		center <- fm.conv.FM2R(ret[[1]])
		all.finite <- fm.conv.FM2R(ret[[2]])
		wx.cs <- fm.conv.FM2R(ret[[3]])
		x.cp <- ret[[4]]
	}
	else {
		ret <- fm.materialize(all.finite, wx.cs, x.cp)
		all.finite <- fm.conv.FM2R(ret[[1]])
		wx.cs <- fm.conv.FM2R(ret[[2]])
		x.cp <- ret[[3]]
	}
	if (!all.finite)
		stop("'x' must contain finite values only")

	x.cp <- x.cp - wx.cs %*% t(center) - center %*% t(wx.cs) + center %*% t(center)
	cov <- switch(match.arg(method), unbiased = x.cp/(1 - sum(wt^2)), ML = x.cp)
	y <- list(cov = cov, center = center, n.obs = n)
	if (with.wt)
		y$wt <- wt
	if (cor) {
		Is <- 1/sqrt(diag(cov))
		R <- cov
		R[] <- Is * cov * rep(Is, each = nrow(cov))
		y$cor <- R
	}
	.set.test.na(orig.test.na)
	y
}

#' Correlation, Variance and Covariance (Matrices)
#'
#' \code{cov} and \code{cor} compute the covariance or correlation of \code{x}
#' and \code{y} if these are vectors.  If \code{x}
#' and \code{y} are matrices then the covariances (or correlations) between
#' the columns of \code{x} and the columns of \code{y} are computed.
#' @param x a numeric vector or matrix
#' @param y \code{NULL} (default) or a vector or matrix with compatible
#' dimensions to \code{x}. The default is equivalent to \code{y=x}.
#' @param use an optional character string giving a method for computing
#'        covariances in the presence of missing values.  This must be
#'        (an abbreviation of) one of the strings \code{"everything"},
#'        \code{"all.obs"}, \code{"complete.obs"}, \code{"na.or.complete"}, or
#'        \code{"pairwise.complete.obs"}.
#' @param method a character string indicating which correlation coefficient
#'        (or covariance) is to be computed.  One of \code{"pearson"}
#'        (default), \code{"kendall"}, or \code{"spearman"}: can be abbreviated.
#'        Right now this argument isn't used.
#' @name cor
NULL

#' @rdname cor
setMethod("cor", "fm", .cor.int)
#' @rdname cor
setMethod("cov", "fm", .cov.int)

#' Weighted Covariance Matrices
#'
#' Returns a list containing estimates of the weighted covariance
#' matrix and the mean of the data, and optionally of the (weighted)
#' correlation matrix.
#'
#' By default, \code{method = "unbiased"}, The covariance matrix is
#' divided by one minus the sum of squares of the weights, so if the
#' weights are the default (1/n) the conventional unbiased estimate
#' of the covariance matrix with divisor (n - 1) is obtained.  This
#' differs from the behaviour in S-PLUS which corresponds to \code{method
#' = "ML"} and does not divide.
#'
#' @param x a matrix. As usual, rows are observations and columns are variables.
#' @param wt a non-negative and non-zero vector of weights for each
#'        observation.  Its length must equal the number of rows of
#'        \code{x}.
#' @param cor a logical indicating whether the estimated correlation
#'        weighted matrix will be returned as well.
#' @param center either a logical or a numeric vector specifying the centers
#'        to be used when computing covariances.  If \code{TRUE}, the
#'        (weighted) mean of each variable is used, if \code{FALSE}, zero is
#'        used.  If \code{center} is numeric, its length must equal the
#'        number of columns of \code{x}.
#' @param method string specifying how the result is scaled, see \code{Details}
#'        below. Can be abbreviated.
#' @return A list containing the following named components:
#' \itemize{
#' \item{cov}{the estimated (weighted) covariance matrix.}
#' \item{center}{an estimate for the center (mean) of the data.}
#' \item{n.obs}{the number of observations (rows) in \code{x}.}
#' \item{wt}{the weights used in the estimation.  Only returned if given
#'          as an argument.}
#' \item{cor}{the estimated correlation matrix.  Only returned if \code{cor} is
#'          \code{TRUE}}
#' }
setMethod("cov.wt", "fm", .cov.wt.int)
