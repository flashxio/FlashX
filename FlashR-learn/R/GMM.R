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

esti.cov.full <- function(X, resp, nk, means, reg.covar)
{
	k <- nrow(means)
	d <- ncol(means)
	covs <- list()
	for (i in 1:k) {
		diff <- sweep(X, 2, means[i,], "-")
		covs[[i]] <- t(resp[,i] * diff) %*% diff / nk[i]
	}
	lapply(covs, function(x) as.matrix(x) + diag(rep(reg.covar, d)))
}

# In this case, we assume all Gaussian distribution has the same covariance matrix.
esti.cov.tied <- function(X, resp, nk, means, reg.covar)
{
	avg.X2 <- t(X) %*% X
	avg.means2 <- t(means * nk) %*% means
	as.matrix((avg.X2 - avg.means2) / sum(nk)) + diag(rep(reg.covar, d))
}

esti.cov.diag <- function(X, resp, nk, means, reg.covar)
{
	avg.X2 <- t(resp) %*% (X * X) / nk
	avg.means2 <- means^2
	avg.Xmeans <- means * (t(resp) %*% X) / nk
	ret <- sweep(as.matrix(avg.X2 - 2*avg.Xmeans + avg.means2), 2, reg.covar, "+")
	ret
}

esti.cov.spherical <- function(X, resp, nk, means, reg.covar)
{
	rowMeans(esti.cov.diag(X, resp, nk, means, reg.covar))
}

# This estimate the parameters of Mixture of Gaussian
esti.gaussian.params <- function(X, resp, reg.covar, cov.type)
{
	print(any(is.nan(resp)))
	n <- nrow(X)
	nk <- colSums(resp)
	# a k x d matrix. Each row is the mean of a component.
	means <- t(resp) %*% X / nk
	# a list of covariance matrices for all components.
	covs <- if (cov.type == "full") esti.cov.full(X, resp, nk, means, reg.covar)
		else if (cov.type == "tied") esti.cov.tied(X, resp, nk, means, reg.covar)
		else if (cov.type == "diag") esti.cov.diag(X, resp, nk, means, reg.covar)
		else if (cov.type == "spherical") esti.cov.spherical(X, resp, nk, means, reg.covar)
		else NULL
	list(weights=as.vector(nk)/n, means=as.matrix(means), covs=covs)
}

init.params <- function(X, k, reg.covar, method, cov.type)
{
	N <- nrow(X)
	if (method == "kmeans") {
		res <- kmeans(X, k)
		resp <- fm.as.sparse.matrix(fm.as.factor(res$cluster))
	}
	else if (method == "random") {
		resp <- fm.runif.matrix(N, k, in.mem=fm.in.mem(X))
		# each row needs to sum up to 1.
		resp <- resp / rowSums(resp)
	}
	else
		stop("unknown init method")

	# Estimate weights, means and covariances
	params <- esti.gaussian.params(X, resp, reg.covar, cov.type)
	list(weights=params$weights, means=params$means, covs=params$covs)
}

est.logprob <- function(X, means, covars, cov.type)
{
	n <- nrow(X)
	d <- ncol(X)
	k <- nrow(means)

	comp.logprob.vec <- function(X, mu, covar.vec) {
		X1 <- sweep(X, 2, mu, "-")
		X2 <- sweep(X1, 2, covar.vec, "/")
		rowSums(X2 * X1)
	}

	if (cov.type == "full") {
		logprob <- list()
		for (i in 1:k)
			logprob[[i]] <- fm.dmvnorm(X, means[i,], covars[[i]], log=TRUE)
		return(do.call(cbind, logprob))
	}
	else if (cov.type == "tied") {
		logprob <- list()
		for (i in 1:k)
			logprob[[i]] <- fm.dmvnorm(X, means[i,], covars, log=TRUE)
		return(do.call(cbind, logprob))
	}
	else if (cov.type == "diag") {
		logprob <- list()
		for (i in 1:k) {
			logprob[[i]] <- comp.logprob.vec(X, means[i,], covars[i,])
			logprob[[i]] <- -sum(log(covars[i,]))/2 - 0.5 * (d * log(2 * pi) + logprob[[i]])
		}
		return(do.call(cbind, logprob))
	}
	else if (cov.type == "spherical") {
		logprob <- list()
		for (i in 1:k) {
			logprob[[i]] <- comp.logprob.vec(X, means[i,], rep(covars[i], d))
			logprob[[i]] <- -log(covars[i])*d/2 - 0.5 * (d * log(2 * pi) + logprob[[i]])
		}
		return(do.call(cbind, logprob))
	}
}

est.weighted.logprob <- function(X, means, covars, cov.type, weights)
{
	sweep(est.logprob(X, means, covars, cov.type), 2, log(weights), "+")
}

logsumexp <- function(X)
{
	max.X <- fm.agg.mat(X, 1, fm.bo.max)
	log(rowSums(exp(X - max.X))) + max.X
}

# Estimate the log likelihood
# @return norm
# @return resp is a n x k matrix. It indicates the probability
#        that a data point belongs to a cluster.
fm.estep <- function(X, params, cov.type)
{
	weighted.logprob <- est.weighted.logprob(X, params$means,
											 params$covs, cov.type,
											 params$weights)
	print(paste("weighted:", any(is.nan(weighted.logprob))))
	logprob.norm <- logsumexp(weighted.logprob)
	log.resp <- weighted.logprob - logprob.norm
	fm.materialize(log.resp, logprob.norm)
	list(norm=as.vector(mean(logprob.norm)), resp=log.resp)
}

# This estimate the parameters.
fm.mstep <- function(X, log.resp, reg.covar, cov.type)
{
	params <- esti.gaussian.params(X, exp(log.resp), reg.covar, cov.type)
	list(weights=params$weights, means=params$means, covs=params$covs)
}

compute.lower.bound <- function(log.resp, log.norm)
{
	log.norm
}

# Fit a mixture of Gaussian distribution on the data.
#
# @param X is a n x d matrix. It's the input data.
# @param k the number of components.
# @param reg.covar is a real value. It's added to the diagonal of
#        the covariance matrix to make it non-singular.
# @param cov.type is the type of covariance matrix. It can be
#        one of {"full", "tied", "diag", "spherical"}.
#        \itemize{
#        \item{"full"}{each component has its own general covariance matrix.}
#        \item{"tied"}{all components share the same general covariance matrix.}
#        \item{"diag"}{each component has its own diagonal covariance matrix.}
#        \item{"spherical"}{each component has its own single variance.}
#        }
# @return 
#        \itemize{
#        \item{loglik}{a n x k matrix, whose \code{[i, k]}th entry is
#                the conditional probability of the ith observation
#                belonging to the kth component of the mixture.}
#        \item{iter}{the number of iterations}
#        \item{parameters}{parameters of the mixture of Gaussian distribution.
#             \itemize{
#             \item{weights}{a vector of k elements. Each element is
#              the weight of the Gaussian distribution in the mixture.}
#             \item{means}{a k x d matrix. Each row is the mean of a Gaussian distribution.}
#             \item{covs}{a list of matrices, a matrix or a vector, depending on \code{cov.type}}}
#        }
GMM.fit <- function(X, k, max.iter=100, tol=1e-3, reg.covar=1e-6,
					method=c("random", "kmeans"),
					cov.type=c("full", "tied", "diag", "spherical"))
{
	method <- match.arg(method)
	cov.type <- match.arg(cov.type)
	params <- init.params(X, k, reg.covar, method, cov.type)
	for (i in 1:max.iter) {
		eret <- fm.estep(X, params, cov.type)
		params <- fm.mstep(X, eret$resp, reg.covar, cov.type)
		lb <- compute.lower.bound(eret$resp, eret$norm)
		print(lb)
		if (i > 5 && abs(lb - prev.lb) < tol)
			break
		prev.lb <- lb
		gc()
	}
	structure(list(loglik=eret$resp, score=eret$norm, iter=i,
				   cov.type=cov.type, parameters=params), class="GMM")
}

get.nparameters <- function(object)
{
	cov.params <- 0
	if (object$cov.type == "full")
		cov.params <- length(object$parameters$covs[[1]]) * length(object$parameters$covs)
	else
		cov.params <- length(object$parameters$covs)
	cov.params + length(object$parameters$means) + length(object$parameters$weights) - 1
}

BIC.GMM <- function(object, ...)
{
	n <- nrow(object$loglik)
	k <- nrow(object$parameters$centers)
	d <- ncol(object$parameters$centers)
	-2 * object$score * n + log(n)*get.nparameters(object)
}

AIC.GMM <- function(object, ..., k=2)
{
	n <- nrow(object$loglik)
	k <- nrow(object$parameters$centers)
	d <- ncol(object$parameters$centers)
	-2 * object$score * n + 2*get.nparameters(object)
}
