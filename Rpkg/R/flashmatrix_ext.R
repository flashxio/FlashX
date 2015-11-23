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
	!fm.mapply2(o1, o2, fm.bo.gt)
}

`<=.fmV`  <- function(o1, o2)
{
	!fm.mapply2(o1, o2, fm.bo.gt)
}

`<.fm`  <- function(o1, o2)
{
	!fm.mapply2(o1, o2, fm.bo.ge)
}

`<.fmV`  <- function(o1, o2)
{
	!fm.mapply2(o1, o2, fm.bo.ge)
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
setMethod("length", signature(x = "fmV"), function(x) x@len)
setMethod("typeof", signature(x = "fm"), function(x) fm.typeof(x))
setMethod("typeof", signature(x = "fmV"), function(x) fm.typeof(x))
setMethod("t", signature(x = "fm"), function(x) fm.t(x))

# Aggregation on a FlashMatrixR vector/matrix.
setMethod("sum", signature(x = "fm"), function(x) fm.agg(x, fm.bo.add))
setMethod("sum", signature(x = "fmV"), function(x) fm.agg(x, fm.bo.add))
setMethod("min", signature(x = "fm"), function(x) fm.agg(x, fm.bo.min))
setMethod("min", signature(x = "fmV"), function(x) fm.agg(x, fm.bo.min))
setMethod("max", signature(x = "fm"), function(x) fm.agg(x, fm.bo.max))
setMethod("max", signature(x = "fmV"), function(x) fm.agg(x, fm.bo.max))

# Miscellaneous Mathematical Functions
setMethod("abs", signature(x = "fm"), function(x) fm.sapply(x, fm.buo.abs))
setMethod("abs", signature(x = "fmV"), function(x) fm.sapply(x, fm.buo.abs))
setMethod("sqrt", signature(x = "fm"), function(x) fm.sapply(x, fm.buo.sqrt))
setMethod("sqrt", signature(x = "fmV"), function(x) fm.sapply(x, fm.buo.sqrt))
setMethod("ceiling", signature(x = "fm"), function(x) fm.sapply(x, fm.buo.ceil))
setMethod("ceiling", signature(x = "fmV"), function(x) fm.sapply(x, fm.buo.ceil))
setMethod("floor", signature(x = "fm"), function(x) fm.sapply(x, fm.buo.floor))
setMethod("floor", signature(x = "fmV"), function(x) fm.sapply(x, fm.buo.floor))

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
