# Copyright 2014 Open Connectome Project (http://openconnecto.me)
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

# This file contains the FlashR functions implemented with
# the FlashR functions in flashmatrix.R.

#' Maxima and Minima
#'
#' Returns the (parallel) maxima and minima of the input values.
#'
#' `fm.pmin2' and `fm.pmax2' works on two input vectors and returns
#' return a single vector/matrix giving the 'parallel' maxima (or minima) of
#' the vectors/matrices. The first element of the result is the maximum (minimum)
#' of the first elements of all the arguments, the second element of
#' the result is the maximum (minimum) of the second elements of all
#' the arguments and so on.
#'
#' @param o1 the input vector/matrix.
#' @param o2 the input vector/matrix.
#' @return a vector/matrix.
#' @name minmax
#' @author Da Zheng <dzheng5@@jhu.edu>
fm.pmin2 <- function(o1, o2)
{
	fm.mapply2(o1, o2, fm.bo.min)
}

#' @name minmax
fm.pmax2 <- function(o1, o2)
{
	fm.mapply2(o1, o2, fm.bo.max)
}

`+.fm` <- function(o1, o2)
{
	fm.mapply2(o1, o2, fm.bo.add)
}

`+.fmV` <- function(o1, o2)
{
	fm.mapply2(o1, o2, fm.bo.add)
}

`-.fm` <- function(o1, o2)
{
	fm.mapply2(o1, o2, fm.bo.sub)
}

`-.fmV` <- function(o1, o2)
{
	fm.mapply2(o1, o2, fm.bo.sub)
}

`*.fm` <- function(o1, o2)
{
	fm.mapply2(o1, o2, fm.bo.mul)
}

`*.fmV` <- function(o1, o2)
{
	fm.mapply2(o1, o2, fm.bo.mul)
}

`/.fm` <- function(o1, o2)
{
	fm.mapply2(o1, o2, fm.bo.div)
}

`/.fmV` <- function(o1, o2)
{
	fm.mapply2(o1, o2, fm.bo.div)
}

`==.fm`  <- function(o1, o2)
{
	fm.mapply2(o1, o2, fm.bo.eq)
}

`==.fmV`  <- function(o1, o2)
{
	fm.mapply2(o1, o2, fm.bo.eq)
}

`!=.fm`  <- function(o1, o2)
{
	!fm.mapply2(o1, o2, fm.bo.eq)
}

`!=.fmV`  <- function(o1, o2)
{
	!fm.mapply2(o1, o2, fm.bo.eq)
}

`>.fm`  <- function(o1, o2)
{
	fm.mapply2(o1, o2, fm.bo.gt)
}

`>.fmV`  <- function(o1, o2)
{
	fm.mapply2(o1, o2, fm.bo.gt)
}

`>=.fm`  <- function(o1, o2)
{
	fm.mapply2(o1, o2, fm.bo.ge)
}

`>=.fmV`  <- function(o1, o2)
{
	fm.mapply2(o1, o2, fm.bo.ge)
}

`<=.fm`  <- function(o1, o2)
{
	fm.mapply2(o1, o2, fm.bo.le)
}

`<=.fmV`  <- function(o1, o2)
{
	fm.mapply2(o1, o2, fm.bo.le)
}

`<.fm`  <- function(o1, o2)
{
	fm.mapply2(o1, o2, fm.bo.lt)
}

`<.fmV`  <- function(o1, o2)
{
	fm.mapply2(o1, o2, fm.bo.lt)
}

`!.fm` <- function(o)
{
	fm.sapply(o, fm.buo.not)
}

`!.fmV` <- function(o)
{
	fm.sapply(o, fm.buo.not)
}

setMethod("%*%", signature(x = "fm", y = "fm"), function(x, y) fm.multiply(x, y))
setMethod("%*%", signature(x = "fm", y = "fmV"), function(x, y) fm.multiply(x, y))

setMethod("dim", signature(x = "fm"), function(x) c(x@nrow, x@ncol))
setMethod("nrow", signature(x = "fm"), function(x) x@nrow)
setMethod("ncol", signature(x = "fm"), function(x) x@ncol)
setMethod("length", signature(x = "fmV"), function(x) x@len)
setMethod("length", signature(x = "fm"), function(x) x@nrow * x@ncol)
setMethod("typeof", signature(x = "fm"), function(x) fm.typeof(x))
setMethod("typeof", signature(x = "fmV"), function(x) fm.typeof(x))
setMethod("t", signature(x = "fm"), function(x) fm.t(x))

# Aggregation on a FlashMatrixR vector/matrix.
setMethod("sum", signature(x = "fm"), function(x) fm.agg(x, fm.bo.add))
setMethod("sum", signature(x = "fmV"), function(x) fm.agg(x, fm.bo.add))
setMethod("mean", signature(x = "fm"), function(x) sum(x)/length(x))
setMethod("mean", signature(x = "fmV"), function(x) sum(x)/length(x))
setMethod("min", signature(x = "fm"), function(x) fm.agg(x, fm.bo.min))
setMethod("min", signature(x = "fmV"), function(x) fm.agg(x, fm.bo.min))
setMethod("max", signature(x = "fm"), function(x) fm.agg(x, fm.bo.max))
setMethod("max", signature(x = "fmV"), function(x) fm.agg(x, fm.bo.max))
setMethod("sd", signature(x = "fmV"), function(x) {
		  n <- length(x)
		  avg <- mean(x)
		  sqrt((sum(x * x) - n * avg * avg) / (n - 1))
})
setMethod("sd", signature(x = "fm"), function(x) {
		  n <- length(x)
		  avg <- mean(x)
		  sqrt((sum(x * x) - n * avg * avg) / (n - 1))
})

# Miscellaneous Mathematical Functions
setMethod("abs", signature(x = "fm"), function(x) fm.sapply(x, fm.buo.abs))
setMethod("abs", signature(x = "fmV"), function(x) fm.sapply(x, fm.buo.abs))
setMethod("sqrt", signature(x = "fm"), function(x) fm.sapply(x, fm.buo.sqrt))
setMethod("sqrt", signature(x = "fmV"), function(x) fm.sapply(x, fm.buo.sqrt))
setMethod("ceiling", signature(x = "fm"), function(x) fm.sapply(x, fm.buo.ceil))
setMethod("ceiling", signature(x = "fmV"), function(x) fm.sapply(x, fm.buo.ceil))
setMethod("floor", signature(x = "fm"), function(x) fm.sapply(x, fm.buo.floor))
setMethod("floor", signature(x = "fmV"), function(x) fm.sapply(x, fm.buo.floor))
setMethod("round", signature(x = "fm"), function(x) fm.sapply(x, fm.buo.round))
setMethod("round", signature(x = "fmV"), function(x) fm.sapply(x, fm.buo.round))
setMethod("log", signature(x = "fm"), function(x, base) fm.sapply(x, fm.buo.log))
setMethod("log", signature(x = "fmV"), function(x, base) fm.sapply(x, fm.buo.log))
setMethod("log10", signature(x = "fm"), function(x) fm.sapply(x, fm.buo.log10))
setMethod("log10", signature(x = "fmV"), function(x) fm.sapply(x, fm.buo.log10))
setMethod("log2", signature(x = "fm"), function(x) fm.sapply(x, fm.buo.log2))
setMethod("log2", signature(x = "fmV"), function(x) fm.sapply(x, fm.buo.log2))
setMethod("exp", signature(x = "fm"), function(x) fm.mapply2(exp(1), x, fm.bo.pow))
setMethod("exp", signature(x = "fmV"), function(x) fm.mapply2(exp(1), x, fm.bo.pow))

setMethod("rowSums", signature(x = "fm", na.rm = "ANY", dims = "ANY"),
		  function(x, na.rm, dims) {
			  fm.agg.mat(x, 1, fm.bo.add)
		  })
setMethod("colSums", signature(x = "fm", na.rm = "ANY", dims = "ANY"),
		  function(x, na.rm, dims) {
			  fm.agg.mat(x, 2, fm.bo.add)
		  })
setMethod("rowMeans", signature(x = "fm", na.rm = "ANY", dims = "ANY"),
		  function(x, na.rm, dims) {
			  rowSums(x, na.rm, dims) / dim(x)[2]
		  })

setMethod("colMeans", signature(x = "fm", na.rm = "ANY", dims = "ANY"),
		  function(x, na.rm, dims) {
			  colSums(x, na.rm, dims) / dim(x)[1]
		  })

# Print FlashMatrix objects.
setMethod("print", signature(x = "fm"), function(x)
		  cat("FlashMatrixR matrix ", x@name, ": ", dim(x)[1], " rows, ", dim(x)[2],
			  " columns, is sparse: ", fm.is.sparse(x), "\n", sep=""))
setMethod("print", signature(x = "fmV"), function(x)
	cat("FlashVectorR vector ", x@name, ": length: ", length(x), "\n", sep=""))
setMethod("print", signature(x = "fm.bo"), function(x)
	cat("FLashMatrixR basic operator:", x@name, "\n"))

fm.table <- function(x)
{
	count <- fm.create.agg.op(fm.bo.count, fm.bo.add, "count")
	fm.sgroupby(x, count)
}

fm.as.integer <- function(x) fm.sapply(x, fm.buo.as.int)
fm.as.numeric <- function(x) fm.sapply(x, fm.buo.as.numeric)
