# Copyright 2014 Open Connectome Project (http://openconnecto.me)
# Written by Da Zheng (zhengda1936@gmail.com)
#
# This file is part of FlashGraph.
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

# This implements Implicit Restart Lanczos Method.

initr <- function(n)
{
	x <- rep.int(1, n)
#	x <- runif(n)
	x <- x / norm(x, type="2")
}

lanczos <- function(A.multiply, k, m, r, V, T)
{
	n <- length(r)
	b <- norm(r, type="2")
	for (j in (k+1):m) {
		restart <- FALSE
		# Step 1: we need to reinitialize r if necessary.
		if (b == 0) {
			r <- initr(n)
			restart <- TRUE
		}

		# Step 2
#		stopifnot(b >= safmin)
		vj <- r / b
		V[,j] <- vj

		# Step 3
		r <- (A.multiply(vj))[,1]

		# Step 4
		wnorm <- norm(r, type="2")
		w <- t(V[,1:j]) %*% r
		r <- r - V[,1:j] %*% w
		a <- w[j]
		T[j,j] <- a
		if (j > 1) {
			if (restart) {
				T[j, j-1] <- 0
				T[j-1, j] <- 0
			}
			else {
				T[j, j-1] <- b
				T[j-1, j] <- b
			}
		}
		b <- norm(r, type="2")

		# Step 5
		orth.iter <- 0
		while (b <= 0.717 * wnorm && orth.iter < 2) {
			orth.iter <- orth.iter + 1
			s <- (t(V[,1:j]) %*% r)
			r <- r - (V[,1:j] %*% s)[,1]
			a <- a + s[j,1]
			T[j,j] <- a
			new.b <- norm(r, type="2")
			if (new.b > 0.717 * b) {
				b <- new.b
				break
			}
			else {
				b <- new.b
				next
			}
			r <- 0
			b <- 0
		}

		if (j > 1 && T[j, j-1] < 0) {
			T[j, j-1] <- -T[j, j-1]
			T[j-1, j] <- -T[j-1, j]
			if (j < m)
				V[,j] <- -V[,j]
			else
				r <- -r
		}
	}
	list(T, V, r, b)
}

compute.eigen <- function(T, which, last.beta, k, m, tol)
{
	ev <- eigen(T, TRUE)
	if (which == "LA") {
		ret <- sort(ev$values, decreasing=TRUE, index.return=TRUE)
	}
	else if (which == "SA") {
		ret <- sort(ev$values, decreasing=FALSE, index.return=TRUE)
	}
	else if (which == "LM") {
		ret <- sort(abs(ev$values), decreasing=TRUE, index.return=TRUE)
	}
	else if (which == "SM") {
		ret <- sort(abs(ev$values), decreasing=FALSE, index.return=TRUE)
	}
	sorted.eigen.values <- ev$values[ret$ix]
	sorted.eigen.vecs <- ev$vectors[,ret$ix]

	m <- dim(sorted.eigen.vecs)[1]
	converged <- c()
	for (i in 1:k) {
		bound <- abs(last.beta * sorted.eigen.vecs[m, i])
		if (bound < tol * abs(sorted.eigen.values[i])) {
			converged <- append(converged, i)
		}
	}

	# wanted eigen values
	list(wval=sorted.eigen.values[converged],
		 # wanted eigen vectors
		 wvec=sorted.eigen.vecs[, converged],
		 # unwanted eigen values
		 uval=sorted.eigen.values[(k+1):length(sorted.eigen.values)])
}

IRLM <- function(A.multiply, n, k, m, which, tol)
{
	p <- m - k
#	n <- dim(A)[1]
	V <- matrix(rep.int(0, n*m), nrow=n, ncol=m)
	T <- matrix(rep.int(0, m*m), nrow=m, ncol=m)
	r <- initr(n)
	res <- lanczos(A.multiply, 0, m, r, V, T)
	T <- res[[1]]
	V <- res[[2]]
	rm <- res[[3]]
	last.beta <- res[[4]]
	while (TRUE) {
		ev <- compute.eigen(T, which, last.beta, k, m, tol)
		if (length(ev$wval) >= k) {
			return(list(values=ev$wval[1:k], vectors=V %*% ev$wvec[, 1:k]))
		}

		Q <- diag(m)
		for (j in 1:p) {
			mu <- ev$uval[j]
			qr.res <- qr(T - mu * diag(m), tol=tol)
			Qj <- qr.Q(qr.res)
			T <- t(Qj) %*% T
			T <- T %*% Qj
			Q <- Q %*% Qj
		}
		bk <- T[k+1, k]
		sigma <- Q[m, k]
		rk <-  rm * sigma
		V[, 1:k] <- V %*% Q[,1:k];
		tmp <- matrix(rep.int(0, m * m), nrow=m, ncol=m)
		tmp[1:k,1:k] <- T[1:k,1:k]
		T <- tmp

		res <- lanczos(A.multiply, k, m, rk, V, T)
		T <- res[[1]]
		V <- res[[2]]
		rm <- res[[3]]
		last.beta <- res[[4]]
	}
}

SVD <- function(A, k, m)
{
	A.multiply <- function(x) {
		t(A) %*% (A %*% x)
	}
	ev <- IRLM(A.multiply, dim(A)[2], k, m)
	sqrt(ev)
}

#IRLM(A.multiply, dim(A)[1], 5, 10)
#SVD(A, 5, 10)
