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
	x <- as.double(x)
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
	sum.x <- sum(x)
	sum.x2 <- sum(x2)

	if (test.na) {
		x.is.na <- any(.is.na.only(x))
		if (.fmV2scalar(x.is.na))
			return(.get.na(typeof(x)))
	}
	else {
		# If we remove NA, we should calculate the number of
		# NAs in the vector.
		sum.na <- sum(in.is.na)
		n <- n - .fmV2scalar(sum.na)
	}
	sum.x <- .fmV2scalar(sum.x)
	sum.x2 <- .fmV2scalar(sum.x2)
	avg <- sum.x / n
	sqrt((sum.x2 - n * avg * avg) / (n - 1))
}

#' Standard Deviation
#'
#' This function computes the standard deviation of the values in \code{x}.
#' If \code{na.rm} is \code{TRUE} then missing values are removed before
#' computation proceeds.
#'
#' @param x a FlashR vector or matrix.
#' @param na.rm logical. Should missing values be removed?
#' @name sd
#'
#' @examples
#' mat <- fm.runif.matrix(100, 10)
#' sd(mat)
NULL

#' @rdname sd
setMethod("sd", "fm", .sd.int)
#' @rdname sd
setMethod("sd", "fmV", .sd.int)

.cov.int <- function(x, y=NULL, use="everything",
				   method=c("pearson", "kendall", "spearman"))
{
	x <- as.double(x)
	if (!is.null(y)) {
		stopifnot(nrow(x) == nrow(y))
		y <- as.double(y)
	}
	n <- nrow(x)
	if (is.null(y)) {
		# we need to transpose x and compute rowSum instead of computing
		# colSum on the original matrix. The reason is that fm.materialize
		# only works on a set of matrices the same the long dimension and
		# dimension size. TODO I need to change that.
		x.sum <- rowSums(t(x))
		x.prod <- fm.multiply(t(x), x)
		x.mu <- x.sum / n
		x.mu <- fm.conv.FM2R(x.mu)
		ret <- (fm.conv.FM2R(x.prod) - n * x.mu %*% t(x.mu)) / (n - 1)
	}
	else {
		x.sum <- rowSums(t(x))
		y.sum <- rowSums(t(y))
		xy.prod <- fm.multiply(t(x), y)
		x.mu <- x.sum / n
		y.mu <- y.sum / n
		x.mu <- fm.conv.FM2R(x.mu)
		y.mu <- fm.conv.FM2R(y.mu)
		ret <- (fm.conv.FM2R(xy.prod) - n * x.mu %*% t(y.mu)) / (n - 1)
	}
	ret
}

.cor.int <- function(x, y=NULL, use="everything",
				   method=c("pearson", "kendall", "spearman"))
{
	x <- as.double(x)
	if (!is.null(y)) {
		stopifnot(nrow(x) == nrow(y))
		y <- as.double(y)
	}
	n <- nrow(x)
	if (is.null(y)) {
		x.sum <- rowSums(t(x))
		x2.sum <- rowSums(t(x * x))
		x.prod <- fm.multiply(t(x), x)
		x.mu <- x.sum / n
		x.sd <- sqrt((x2.sum - n * x.mu * x.mu) / (n - 1))
		x.mu <- fm.conv.FM2R(x.mu)
		x.sd <- fm.conv.FM2R(x.sd)
		ret <- (fm.conv.FM2R(x.prod) - n * x.mu %*% t(x.mu)) / (n - 1) / (x.sd %*% t(x.sd))
	}
	else {
		x.sum <- rowSums(t(x))
		y.sum <- rowSums(t(y))
		x2.sum <- rowSums(t(x * x))
		y2.sum <- rowSums(t(y * y))
		xy.prod <- fm.multiply(t(x), y)
		x.mu <- x.sum / n
		x.sd <- sqrt((x2.sum - n * x.mu * x.mu) / (n - 1))
		x.mu <- fm.conv.FM2R(x.mu)
		x.sd <- fm.conv.FM2R(x.sd)
		y.mu <- y.sum / n
		y.sd <- sqrt((y2.sum - n * y.mu * y.mu) / (n - 1))
		y.mu <- fm.conv.FM2R(y.mu)
		y.sd <- fm.conv.FM2R(y.sd)
		ret <- (fm.conv.FM2R(xy.prod) - n * x.mu %*% t(y.mu)) / (n - 1) / (x.sd %*% t(y.sd))
	}
	ret
}

.cov.wt.int <- function (x, wt = rep(1/nrow(x), nrow(x)), cor = FALSE, center = TRUE,
					       method = c("unbiased", "ML"))
{
	if (is.data.frame(x))
		x <- fm.as.matrix(x)
	else if (!fm.is.matrix(x))
		stop("'x' must be a matrix or a data frame")
	x <- as.double(x)
	n <- nrow(x)
	if (with.wt <- !missing(wt)) {
		if (length(wt) != n)
			stop("length of 'wt' must equal the number of rows in 'x'")
		any.neg <- fm.any(wt < 0, TRUE)
		s <- sum(wt)
		if (any.neg[1] || s[1] == 0)
			stop("weights must be non-negative and not all zero")
		wt <- wt/s
	}

	all.finite <- fm.all(is.finite(t(x)), TRUE)
	if (is.logical(center)) {
		center <- if (center)
			rowSums(t(wt * x))
		else 0
	}
	else {
		if (length(center) != ncol(x))
			stop("length of 'center' must equal the number of columns in 'x'")
	}
	wx <- wt * x
	wx.cs <- rowSums(t(wx))
	x.cp <- crossprod(wx, x)
	if (fm.is.sink(center)) {
		center <- fm.conv.FM2R(center)
		all.finite <- fm.conv.FM2R(all.finite)
		wx.cs <- fm.conv.FM2R(wx.cs)
		x.cp <- x.cp
	}
	else {
		all.finite <- fm.conv.FM2R(all.finite)
		wx.cs <- fm.conv.FM2R(wx.cs)
		x.cp <- x.cp
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
	y
}

#' Correlation, Variance and Covariance (Matrices)
#'
#' \code{cov} and \code{cor} compute the covariance or correlation of \code{x}
#' and \code{y} if these are vectors.  If \code{x}
#' and \code{y} are matrices then the covariances (or correlations) between
#' the columns of \code{x} and the columns of \code{y} are computed.
#'
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
#'
#' @examples
#' mat <- fm.runif.matrix(100, 10)
#' var(mat)
#' cor(mat)
#' cov(mat)
NULL

#' @rdname cor
setMethod("var", "fm", function(x, y=NULL, na.rm=FALSE, use) {
		  if (missing(use))
			  use <- if (na.rm) "na.or.complete" else "everything"
		  cor(x, y, use)
})
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
#' @name cov.wt
#'
#' @examples
#' mat <- fm.runif.matrix(100, 10)
#' cov.wt(mat, fm.runif(nrow(mat)))
setMethod("cov.wt", "fm", .cov.wt.int)
