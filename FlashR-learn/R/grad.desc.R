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

gradient.descent <- function(X, y, get.grad, get.hessian, cost, params)
{
	# Add x_0 = 1 as the first column
	m <- if(is.vector(X)) length(X) else nrow(X)
	if(is.vector(X) || (!all(X[,1] == 1))) X <- cbind(fm.as.matrix(fm.rep.int(1, m)), X)

	num.features <- ncol(X)

	# Initialize the parameters
	theta <- matrix(rep(0, num.features), nrow=1)

	L2 <- function(X) {
		sqrt(sum(X * X))
	}

	method <- params$method
	if (method == "RNS")
		rnorms <- rowSums(X * X)

	eta.est <- 1
	# Look at the values over each iteration
	theta.path <- theta
	for (i in 1:params$num.iters) {
		print(paste("iter:", i))
		g <- get.grad(X, y, theta)
		l <- cost(X, y, theta)
		if (!is.null(get.hessian) && i > 5) {
			n <- nrow(X)
			s <- params$hessian_size * n
			if (method == "Newton") {
				D2 <- get.hessian(X, y, theta)
				H <- t(X) %*% sweep(X, 1, fm.as.vector(D2), FUN="*")
			}
			else if (method == "LS") {
			}
			else if (method == "RNS") {
				D2 <- get.hessian(X, y, theta)
				p <- D2 * rnorms
				p <- p / sum(p)
				# TODO
				q <- pmin(fm.rep.int(1, length(p)), p * s)
				# TODO we need to fix this.
				# If idx is empty, it gets some weird error.
				idx <- which(fm.conv.FM2R(fm.runif(n) < q))
				q <- fm.materialize(q)
				p.sub <- q[idx]
				X.sub <- X[idx, ]
				# TODO
				D2 <- fm.materialize(D2)
				D2.sub <- D2[idx]
				# TODO
				D2.sub <- fm.materialize(D2.sub)
				H <- t(X.sub) %*% sweep(X.sub, 1, fm.as.vector(D2.sub / p.sub), FUN="*")
			}
			else if (method == "Uniform") {
				# TODO we need to fix this.
				idx <- which(fm.conv.FM2R(fm.runif(n) < s/n))
				X.sub <- X[idx, ]
				D2.sub <- get.hessian(X.sub, y[idx], theta)
				H <- t(X.sub) %*% sweep(X.sub, 1, fm.as.vector(D2.sub * n), FUN="*") / s
			}
			H <- as.matrix(H)
			h.max <- max(abs(H[H != 0]))
			if (h.max > 1)
				H <- H / h.max
			g <- as.matrix(g / length(y))
			z <- pcg(H, as.vector(-g), maxiter=1000, tol=1e-06)
			z <- as.matrix(t(z))
			params$linesearch <- TRUE
		}
		else {
			params$linesearch <- TRUE
			g <- as.matrix(g / length(y))
			z <- -g
		}
		l <- as.vector(l)/length(y)
		cat(i,  ": L2(g) =", L2(g), ", cost:", l, "\n")

		eta <- eta.est
		if (params$linesearch) {
			delta = params$c * z %*% t(g)
			while (TRUE) {
				costs <- list()
				costs[[1]] <- cost(X, y, theta + eta * z)
				costs[[2]] <- cost(X, y, theta + eta * params$ro * z)
				costs <- lapply(fm.materialize.list(costs), function(o) as.vector(o))
				if (as.vector(costs[[1]]) / length(y) >= l + delta * eta) eta <- eta * params$ro else break
				if (as.vector(costs[[2]]) / length(y) >= l + delta * eta) eta <- eta * params$ro else break
			}
			print(eta)
		}
		eta.est <- eta * 10
		if (eta == 1)
			params$linesearch <- FALSE

		theta <- theta + z * eta
		if(all(is.na(theta))) break
		theta.path <- rbind(theta.path, theta)
	}

	if(params$out.path) return(theta.path) else return(theta.path[nrow(theta.path),])
}
