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

lanczos <- function(A.multiply, k, m, r, V, T)
{
	b <- norm(r, type="2")
	if (k > 0) {
		T[k, k+1] <- b
		T[k+1, k] <- b
	}
	for (j in (k+1):m) {
		vj <- r / b
		V[,j] <- vj
		r <- (A.multiply(vj))[,1]
		vj_1 <- V[, j - 1]
		a <- (t(r) %*% vj)[1,1]
		r <- r - vj * a
		if (j > 1) {
			r <- r - V[, j - 1] * b
		}
		new.b <- norm(r, type="2")
		if (new.b < sqrt(a * a + b * b)) {
			s <- (t(V[,1:j]) %*% r)
			r <- r - (V[,1:j] %*% s)[,1]
			a <- a + s[j,1]
			new.b <- norm(r, type="2")
		}
		b <- new.b
		T[j,j] <- a
		if (j < m) {
			T[j, j+1] <- b
			T[j+1, j] <- b
		}
	}
	list(T, V, r)
}

IRLM <- function(A.multiply, n, k, m)
 {
	p <- m - k
#	n <- dim(A)[1]
	V <- matrix(rep.int(0, n*m), nrow=n, ncol=m)
	T <- matrix(rep.int(0, m*m), nrow=m, ncol=m)
	r <- rep.int(1, n)
#	r <- runif(n)
	res <- lanczos(A.multiply, 0, m, r, V, T)
	T <- res[[1]]
	V <- res[[2]]
	rm <- res[[3]]
	for (i in 1:5) {
		ev <- eigen(T, TRUE, only.values=TRUE)
		eigen.values <- sort(ev$values, decreasing=TRUE)
		Q <- diag(m)
		for (j in 1:p) {
			mu <- eigen.values[j + k]
			qr.res <- qr(T - mu * diag(m))
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
	}
	eigen.values[1:k]
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
