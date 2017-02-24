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

#' Compute PageRank on a graph.
#'
#' This function assumes a graph is stored in a sparse matrix and computes
#' exact PageRank with a sequence of sparse matrix multiplication.
#'
#' @param A a sparse matrix that stores the graph.
#' @param d a damping factor.
#' @param max.niters the maximal number of iterations.
#' @param epsilon determines the precision of PageRank values.
#' @param verbose indicates whether to print extra information.
#' @return a vector that contains the PageRank value for each vertex.
#' @author Da Zheng <dzheng5@@jhu.edu>
fm.PageRank <- function(A, d = 0.15, max.niters = 30, epsilon = 1e-2,
					 verbose = FALSE)
{
	fm.set.test.na(FALSE)
	N <- dim(A)[1]
	epsilon <- epsilon / N
	cat("There are", N, "vertices.", "\n")
	one <- fm.rep.int(1, N)
	out.deg <- A %*% one
	one <- NULL
	pr1 <- fm.rep.int(1/N, N)
	converge <- 0

	niters <- 0
	A <- t(A)
	while (converge < N && niters < max.niters) {
		start <- Sys.time()
		pr2 <- d/N+(1-d)*(A %*% (pr1/out.deg))
		gc()
		converge <- as.vector(sum(as.double(abs(pr1-pr2) < epsilon)))
		if (verbose)
			cat("iter", niters, ", #converge:", converge, "\n")
		pr1 <- pr2
		niters <- niters + 1
		end <- Sys.time()
		print(end - start)
	}
#	fm.set.test.na(orig.test.na)
	pr2
}
