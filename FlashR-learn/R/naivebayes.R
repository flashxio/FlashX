# Copyright 2017 Open Connectome Project (http://openconnecto.me)
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

naiveBayes.train <- function(x, y)
{
	# For the training process, we need to estimate the mean and
	# standard deviation of each feature.
	y <- fm.as.factor(y)
	fm.set.test.na(FALSE)
	gsum <- fm.groupby(x, 2, y, fm.bo.add)
	gsum2 <- fm.groupby(x * x, 2, y, fm.bo.add)
	fm.materialize(gsum, gsum2)
	cnt <- fm.table(y)
	gmean <- gsum / cnt@Freq
	gsd <- sqrt((gsum2 - cnt@Freq * gmean ^ 2) / (cnt@Freq - 1))
	fm.set.test.na(TRUE)
	list(apriori=as.vector(cnt@Freq), mean=as.matrix(gmean),
		 sd=as.matrix(gsd))
}

naiveBayes.predict <- function(object, newdata)
{
	fm.set.test.na(FALSE)
	apriori <- object$apriori
	mean <- object$mean
	sd <- object$sd
	ncls <- length(apriori)
	loglikely <- list()
	for (i in 1:ncls) {
		# Here is the compute the log likelihood of the data
		# given the normal distribution (mean and sd).
		x <- sweep(newdata, 2, mean[i,], "-") ^ 2
		x <- sweep(x, 2, 2 * sd[i,]^2, "/")
		loglikely <- c(loglikely, log(apriori[i]) - sum(log(sqrt(2 * pi * sd[i, ]^2))) - rowSums(x))
	}
	loglikely <- fm.cbind.list(loglikely)
	cls <- fm.materialize(fm.agg.mat(loglikely, 1, fm.bo.which.max))
	fm.set.test.na(TRUE)
	cls
}
