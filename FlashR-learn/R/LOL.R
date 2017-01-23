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

# m is a D x I matrix. Each column is a training instance and each instance has D features.
# labels contains I labels, one for each training instance
# k specifies #columns of the output matrix
# The output is a D x k matrix.
LOL <- function(m, labels, k, type=c("svd", "rand_dense", "rand_sparse")) {
	counts <- as.data.frame(table(fm.conv.FM2R(labels)))
	num.labels <- length(counts$Freq)
	num.features <- dim(m)[1]
	nv <- k - (num.labels - 1)
	gr.sum <- fm.groupby(m, 1, fm.as.factor(labels, 2), fm.bo.add)
	gr.mean <- fm.mapply.row(gr.sum, counts$Freq, fm.bo.div, FALSE)
	diff <- fm.get.cols(gr.mean, 1) - fm.get.cols(gr.mean, 2)
	diff.l2 <- sqrt(sum(diff * diff))
	stopifnot(diff.l2 != 0)
	diff <- diff / diff.l2
	if (nv == 0)
		return(diff)
	if (type == "svd") {
		# compute class conditional mean.
		# Here I assume dimension size is larger than the number of samples.
		# TODO I need to improve groupby to subtract the class conditional mean.
		rlabels <- fm.conv.FM2R(labels) + 1
		fm.materialize(gr.mean)
		mean.list <- list()
		for (i in 1:length(rlabels))
			mean.list[[i]] <- gr.mean[,rlabels[i]]
		mean.mat <- fm.cbind.list(mean.list)

		svd <- fm.svd(m - mean.mat, nv=nv, nu=nv)
		cbind(diff, svd$u)
	}
	else if (type == "rand_dense")
		cbind(diff, fm.rnorm.matrix(length(diff), nv))
	else if (type == "rand_sparse")
		cbind(diff, fm.rsparse.proj(length(diff), nv, 1/sqrt(length(diff))))
	else {
		print("wrong type")
		NULL
	}
}

QOQ <- function(m, labels, k)
{
	rlabels <- fm.conv.FM2R(labels)
	# we have to make sure they are binary labels and their values
	# are 0 or 1.
	stopifnot(sum(rlabels == 0 | rlabels == 1) == length(rlabels))

	counts <- as.data.frame(table(rlabels))
	print(counts)
	num.labels <- length(counts$Freq)
	num.features <- dim(m)[1]
	nv <- k - (num.labels - 1)
	nv0 <- min(nv, counts$Freq[1])
	nv1 <- min(nv, counts$Freq[2])

	gr.sum <- fm.groupby(m, 1, fm.as.factor(labels, 2), fm.bo.add)
	gr.mean <- fm.mapply.row(gr.sum, counts$Freq, fm.bo.div, FALSE)
	diff <- fm.get.cols(gr.mean, 1) - fm.get.cols(gr.mean, 2)
	diff.l2 <- sqrt(sum(diff * diff))
	stopifnot(diff.l2 != 0)
	diff <- diff / diff.l2
	if (nv == 0)
		return(diff)

	c0.idxs <- which(rlabels == 0)
	c1.idxs <- which(rlabels == 1)
	m0 <- sweep(m[, c0.idxs], 1, gr.mean[,1], "-")
	m1 <- sweep(m[, c1.idxs], 1, gr.mean[,2], "-")
	svd0 <- fm.svd(m0, nv=nv0, nu=nv0)
	svd1 <- fm.svd(m1, nv=nv1, nu=nv1)
	vals <- c(svd0$d, svd1$d)
	vecs <- cbind(diff, svd0$u, svd1$u)
	sort.res <- sort(vals, decreasing=TRUE, index.return=TRUE)
	print(sort.res$x)
	print(sort.res$ix)
	col.idxs <- c(1, sort.res$ix + 1)
	vecs[,col.idxs]
}
