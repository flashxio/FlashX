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

for (cl in fm.cls) {
	setMethod("sd", cl, function(x, na.rm) {
			  n <- length(x)
			  x2 <- x * x
			  test.na <- TRUE
			  num.na <- 0
			  if (na.rm) {
				  zero <- get.zero(typeof(x))
				  in.is.na <- is.na(x)
				  x <- ifelse(in.is.na, zero, x)
				  x2 <- ifelse(in.is.na, zero, x2)
				  test.na <- FALSE
			  }
			  sum.x <- fm.agg.lazy(x, fm.bo.add)
			  sum.x2 <- fm.agg.lazy(x2, fm.bo.add)

			  if (test.na) {
				  x.is.na <- fm.agg.lazy(fm.is.na.only(x), fm.bo.or)
				  res <- fm.materialize(sum.x, sum.x2, x.is.na)
				  if (fmV2scalar(res[[3]]))
					  return(get.na(typeof(x)))
				  sums <- res[1:2]
			  }
			  else {
				  # If we remove NA, we should calculate the number of
				  # NAs in the vector.
				  sum.na <- fm.agg.lazy(in.is.na, fm.bo.add)
				  res <- fm.materialize(sum.x, sum.x2, sum.na)
				  n <- n - fmV2scalar(res[[3]])
				  sums <- res[1:2]
			  }
			  sum.x <- fmV2scalar(sums[[1]])
			  sum.x2 <- fmV2scalar(sums[[2]])
			  avg <- sum.x / n
			  sqrt((sum.x2 - n * avg * avg) / (n - 1))
		  })
}

fm.cov <- function(x, y=NULL, use="everything",
				   method=c("pearson", "kendall", "spearman"))
{
	x.mu <- colSums(x) / nrow(x)
	x0 <- fm.mapply.row(x, x.mu, fm.bo.sub)
	if (is.null(y))
		(t(x0) %*% x0) / (nrow(x) - 1)
	else {
		y.mu <- colSums(y) / nrow(y)
		y0 <- fm.mapply.row(y, y.mu, fm.bo.sub)
		(t(x0) %*% y0) / (nrow(x) - 1)
	}
}

fm.cov.wt <- function (x, wt = rep(1/nrow(x), nrow(x)), cor = FALSE, center = TRUE,
					       method = c("unbiased", "ML"))
{
	if (is.data.frame(x))
		x <- as.matrix(x)
	else if (!is.matrix(x))
		stop("'x' must be a matrix or a data frame")
	if (!all(is.finite(x)))
		stop("'x' must contain finite values only")
	n <- nrow(x)
	if (with.wt <- !missing(wt)) {
		if (length(wt) != n)
			stop("length of 'wt' must equal the number of rows in 'x'")
		if (any(wt < 0) || (s <- sum(wt)) == 0)
			stop("weights must be non-negative and not all zero")
		wt <- wt/s
	}
	if (is.logical(center)) {
		center <- if (center)
			colSums(wt * x)
		else 0
	}
	else {
		if (length(center) != ncol(x))
			stop("length of 'center' must equal the number of columns in 'x'")
	}
	x <- sqrt(wt) * sweep(x, 2, center, check.margin = FALSE)
	cov <- switch(match.arg(method),
				  unbiased = crossprod(x)/(1 - sum(wt^2)), ML = crossprod(x))
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

setMethod("cov", "fm",  fm.cov)
setMethod("cov.wt", "fm", fm.cov.wt)

fm.dmvnorm <- function(X, mu, covar, log=FALSE)
{
	if (fm.is.matrix(covar))
		covar <- fm.conv.FM2R(covar)
	covar.inv <- solve(covar)
	X1 <- sweep(X, 2, mu, "-")
	X2 <- X1 %*% covar.inv
	X3 <- fm.agg.mat(X2 * X1, 1, "+")
	k <- dim(covar)[1]
	dec <- tryCatch(chol(covar), error = function(e) e)
	logdet <- -sum(log(diag(dec)))
	logret <- logdet - 0.5 * k * log(2 * pi) - 0.5 * X3
	if (log)
		ret <- logret
	else
		ret <- exp(logret)
}

fm.summary <- function(x)
{
	lazy.res <- list()
	if (is.matrix(x)) {
		lazy.res[[1]] <- fm.agg.mat.lazy(x, 2, fm.bo.min)
		lazy.res[[2]] <- fm.agg.mat.lazy(x, 2, fm.bo.max)
		lazy.res[[3]] <- fm.agg.mat.lazy(x, 2, fm.bo.add)
		lazy.res[[4]] <- fm.agg.mat.lazy(abs(x), 2, fm.bo.add)
		lazy.res[[5]] <- fm.agg.mat.lazy(x * x, 2, fm.bo.add)
		lazy.res[[6]] <- fm.agg.mat.lazy(x != 0, 2, fm.bo.add)
	}
	else {
		lazy.res[[1]] <- fm.agg.lazy(x, fm.bo.min)
		lazy.res[[2]] <- fm.agg.lazy(x, fm.bo.max)
		lazy.res[[3]] <- fm.agg.lazy(x, fm.bo.add)
		lazy.res[[4]] <- fm.agg.lazy(abs(x), fm.bo.add)
		lazy.res[[5]] <- fm.agg.lazy(x * x, fm.bo.add)
		lazy.res[[6]] <- fm.agg.lazy(x != 0, fm.bo.add)
	}
	res <- fm.materialize.list(lazy.res)
	res <- lapply(res, function(o) fm.conv.FM2R(o))
	mean <- res[[3]]/nrow(x)
	var <- (res[[5]]/nrow(x) - mean^2) * nrow(x) / (nrow(x) - 1)
	list(min=res[[1]], max=res[[2]], mean=mean, normL1=res[[4]],
		 normL2=sqrt(res[[5]]), numNonzeros=res[[6]], var=var)
}
