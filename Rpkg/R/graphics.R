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

fm.hist <- function(x, breaks=c("even", "exp"), exp.base=2, plot=TRUE,
					plot.file="", log="", xlab=NULL, ylab=NULL)
{
	breaks <- match.arg(breaks)
	if (breaks == "exp") {
		x <- ifelse(x > 0, floor(log(x, exp.base)), 0)
		r <- range(x)
		vals <- seq(r[1], r[2], 1)
	}
	else {
		r <- range(x)
		size <- (r[2] - r[1]) / 100
		x <- floor(x / size)
	}
	tbl <- fm.table(x)

	if (breaks == "exp") {
		vals <- fm.conv.FM2R(tbl@val)
		breaks <- exp.base^vals
		x.names <- paste(as.character(exp.base), "^",as.character(vals), sep="")
	}
	else {
		breaks <- fm.conv.FM2R(tbl@val) * size
		x.names <- as.character(breaks)
	}

	if (plot && plot.file != "") {
		pdf(plot.file)
		barplot(fm.conv.FM2R(tbl@Freq), space=0, names.arg=x.names, log=log,
				xlab=xlab, ylab=ylab)
		dev.off()
	}
	else if (plot) {
		barplot(fm.conv.FM2R(tbl@Freq), space=0, names.arg=x.names, log=log,
				xlab=xlab, ylab=ylab)
	}
	else
		list(breaks=breaks, counts=fm.conv.FM2R(tbl@Freq))
}

fm.heatmap <- function(coo1, num.slots, logval=FALSE, plot=TRUE, plot.file="",
					   plot.format="pdf")
{
	tbl1 <- fm.table(floor(coo1[,1] * num.slots))
	tbl2 <- fm.table(floor(coo1[,2] * num.slots))
	min.val1 <- as.vector(min(tbl1@val))
	min.val2 <- as.vector(min(tbl2@val))
	len1 <- as.vector(tbl1@val[length(tbl1@val)]) - as.vector(tbl1@val[1]) + 1
	len2 <- as.vector(tbl2@val[length(tbl2@val)]) - as.vector(tbl2@val[1]) + 1
	vec <- rep(0, times=len1 * len2)
	tbl <- fm.table(floor(coo1[,1] * num.slots)*len2 + floor(coo1[,2]*num.slots))
	min.val <- min.val1 * len2 + min.val2
	vec[as.vector(tbl@val) - min.val + 1] <- as.vector(tbl@Freq)
	row.names <- as.character(seq(min.val2, as.vector(max(tbl2@val)), 1)/num.slots)
	col.names <- as.character(seq(min.val1, as.vector(max(tbl1@val)), 1)/num.slots)
	data <- matrix(vec, len2, len1, dimnames=list(row.names, col.names))
	if (plot) {
		if (plot.file != "" && plot.format == "pdf")
			pdf(plot.file)
		else if (plot.file != "" && plot.format == "jpg")
			jpeg(plot.file)
		if (logval)
			data <- ifelse(log(data) == -Inf, 0, log(data))
		heatmap.2(data, Colv="NA", Rowv="NA", trace="none")
		if (plot.file != "")
			dev.off()
	}
	else {
		data
	}
}
