# Copyright 2014 Open Connectome Project (http://openconnecto.me)
# Written by Da Zheng (zhengda1936@gmail.com)
#
# This file is part of FlashMatrixR.
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

safe.min <- 1e-12

norm2 <- function(x)
{
	sqrt(fm.sum(x * x))
}

initr <- function(n)
{
#	x <- rep.int(1, n)
	x <- fm.runif(n, -1, 1)
}

# In this version of initialization of r, we need to reorthoganize r against
# the current Lanczos basis.
reinitr <- function(n, V)
{
	r <- initr(n)
	b0 <- norm2(r)

	iter <- 0
	while (iter < 5) {
		iter <- iter + 1
		s <- fm.multiply(fm.t(V), r)
		r <- r - fm.multiply(V, s)
		b <- norm2(r)
		if (b <= 0.717 * b0)
			break
		b0 <- b
		print("reothoganize r in reinitr")
	}
	stopifnot(iter < 5)
	r
}

lanczos <- function(A.multiply, k, m, r, V, T)
{
	n <- fm.length(r)
	b <- norm2(r)
	for (j in (k+1):m) {
		restart <- FALSE
		# Step 1: we need to reinitialize r if necessary.
		if (b == 0) {
			stopifnot(j > 1)
			print("reinitialize r")
			r <- reinitr(n, V[,1:(j-1)])
			restart <- TRUE
		}

		# Step 2
#		stopifnot(b >= safmin)
		vj <- r / b
		fm.set.cols(V, j, vj)

		# Step 3
		r <- A.multiply(vj)

		# Step 4
		wnorm <- norm2(r)
		w <- fm.multiply(fm.t(fm.get.cols(V,1:j)), r)
		r <- r - fm.multiply(fm.get.cols(V,1:j), w)
		a <- fm.conv.FM2R(w)[j]
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
		b <- norm2(r)

		# Step 5
		orth.iter <- 0
		while (b <= 0.717 * wnorm && orth.iter < 2) {
			orth.iter <- orth.iter + 1
			s <- fm.multiply(fm.t(fm.get.cols(V,1:j)), r)
			r <- r - fm.multiply(fm.get.cols(V,1:j), s)
			a <- a + fm.conv.FM2R(s)[j]
			T[j,j] <- a
			new.b <- norm2(r)
			if (new.b > 0.717 * b) {
				b <- new.b
				break
			}
			# Try to orthoganize r again
			else if (orth.iter < 2) {
				print("reothoganize r again")
				b <- new.b
				next
			}
			# If we still can't orthoganize r, we need to go to the first
			# step and reinitialize r.
			else {
				print("fail to reorthoganize r")
				r <- 0
				b <- 0
				break
			}
		}

		if (j > 1 && T[j, j-1] < 0) {
			print("negative T[j,j-1]")
			# TODO
#			T[j, j-1] <- -T[j, j-1]
#			T[j-1, j] <- -T[j-1, j]
#			if (j < m)
#				V[,j] <- -V[,j]
#			else
#				r <- -r
		}
	}

	list(T, V, r, b)
}

compute.eigen <- function(T, last.beta)
{
	m <- dim(T)[1]
	ev <- eigen(T, TRUE)
#	print(ev$values)
	bounds <- rep.int(0, m)
	for (i in 1:m)
		bounds[i] <- abs(last.beta * ev$vectors[m, i])
	list(values=ev$values, vectors=ev$vectors, bounds=bounds)
}

sort.eigen <- function(ev, which, k)
{
	if (which == "LA") {
		sort.ret <- sort(ev$values, decreasing=TRUE, index.return=TRUE)
	}
	else if (which == "SA") {
		sort.ret <- sort(ev$values, decreasing=FALSE, index.return=TRUE)
	}
	else if (which == "LM") {
		sort.ret <- sort(abs(ev$values), decreasing=TRUE, index.return=TRUE)
	}
	else if (which == "SM") {
		sort.ret <- sort(abs(ev$values), decreasing=FALSE, index.return=TRUE)
	}
	ret <- list(values=ev$values[sort.ret$ix], vectors=ev$vectors[,sort.ret$ix],
		 bounds=ev$bounds[sort.ret$ix])
	m <- length(ev$values)
	if (m - k > 1) {
		sort.ret <- sort(abs(ret$bounds[(k+1):m]), decreasing=TRUE, index.return=TRUE)
		ret$values[(k+1):m] <- ret$values[k + sort.ret$ix]
		ret$vectors[,(k+1):m] <- ret$vectors[,k + sort.ret$ix]
		ret$bounds[(k+1):m] <- ret$bounds[k + sort.ret$ix]
	}
	ret
}

test.conv <- function(ev, k, tol)
{
	converged <- c()
	for (i in 1:k) {
		if (ev$bounds[i] < tol * abs(ev$values[i])) {
			converged <- append(converged, i)
		}
	}

	m <- length(ev$values)
	# wanted eigen values
	list(wval=ev$values[converged],
		 # wanted eigen vectors
		 wvec=ev$vectors[, converged],
		 # unwanted eigen values
		 uval=ev$values[(k+1):m])
}

#' Compute eigenvalues/vectors on a symmetric matrix.
#'
#' This method reimplements implicitly restarted Lanczos method to compute
#' eigenvalues and eigenvectors. The method is implemented completely in R
#' and follows the implementation of ARPACK library.
#'
#' @param A.multiply the function to perform matrix vector multiplication.
#' @param n Numeric scalar. The dimension of the eigenproblem.
#' @param k Numeric scalar. The number of eigenvalues to be computed.
#' @param m Number of Lanczos vectors to be generated.
#' @param which Specify which eigenvalues/vectors to compute, character
#' constant with exactly two characters.
#' Possible values for symmetric input matrices:
#' "LA" Compute "nev" largest (algebraic) eigenvalues.
#' "SA" Compute "nev" smallest (algebraic) eigenvalues.
#' "LM" Compute "nev" largest (in magnitude) eigenvalues.
#' "SM" Compute "nev" smallest (in magnitude) eigenvalues.
#' @param Numeric scalar. Stopping criterion: the relative accuracy of
#' the Ritz value is considered acceptable if its error is less than 'tol'
#' times its estimated value. If this is set to zero then machine precision is used.
#' @author Da Zheng <dzheng5@@jhu.edu>
fm.arpack <- function(A.multiply, n, k, m, which, tol)
{
#	n <- dim(A)[1]
	V <- fm.matrix(fm.rep.int(0, n*m), nrow=n, ncol=m)
	T <- matrix(rep.int(0, m*m), nrow=m, ncol=m)
	r <- initr(n)
	res <- lanczos(A.multiply, 0, m, r, V, T)
	T <- res[[1]]
	V <- res[[2]]
	rm <- res[[3]]
	last.beta <- res[[4]]
	iter <- 0
	while (TRUE) {
		iter <- iter + 1

		ev <- compute.eigen(T, last.beta)
		ev <- sort.eigen(ev, which, k)
		res <- test.conv(ev, k, tol)
		nconv <- length(res$wval)
		bounds <- ev$bounds
		np <- m - k
		nev <- k
		for (i in (k+1):m) {
			if (bounds[i] == 0) {
				np <- np - 1
				nev <- nev + 1
			}
		}
		
		if (nconv >= k || np == 0) {
			vecs <- fm.conv.R2FM(res$wvec[, 1:k])
			return(list(values=res$wval[1:k], vectors=fm.multiply(V, vecs)))
		}

		if (nconv < nev) {
			nevbef <- nev
			nev <- nev + min(nconv, np / 2)
			if (nev == 1 && m >= 6)
				nev <- m / 2
			else if (nev == 1 && m > 2)
				nev <- 2
			np <- m - nev
			if (nevbef < nev) {
				ev <- sort.eigen(ev, which, nev)
				print("sort eigen again")
			}
			print(nev)
		}

		Q <- diag(m)
		for (j in 1:np) {
			mu <- res$uval[j]
			qr.res <- qr(T - mu * diag(m), tol=tol)
			Qj <- qr.Q(qr.res)
			T <- t(Qj) %*% T %*% Qj
			Q <- Q %*% Qj
		}
		bk <- T[nev+1, nev]
		sigma <- Q[m, nev]
		tmp <- fm.as.vector(fm.get.cols(V,nev+1))
		rk <- bk * tmp + rm * sigma
#		rk <- rm * sigma
		fm.set.cols(V, 1:nev, fm.multiply(V, fm.conv.R2FM(Q[,1:nev])))
		tmp <- matrix(rep.int(0, m * m), nrow=m, ncol=m)
		tmp[1:nev,1:nev] <- T[1:nev,1:nev]
		T <- tmp

		res <- lanczos(A.multiply, nev, m, rk, V, T)
		T <- res[[1]]
		V <- res[[2]]
		rm <- res[[3]]
		last.beta <- res[[4]]
	}
}
