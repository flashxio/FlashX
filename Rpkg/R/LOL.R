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
		fm.cbind(diff, svd$u)
	}
	else if (type == "rand_dense")
		fm.cbind(diff, fm.rnorm.matrix(length(diff), nv))
	else if (type == "rand_sparse")
		fm.cbind(diff, fm.rsparse.proj(length(diff), nv, 1/sqrt(length(diff))))
	else {
		print("wrong type")
		NULL
	}
}
