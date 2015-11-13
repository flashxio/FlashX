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
	fm.ele.wise.op(fm.bo.min, o1, o2)
}

#' @name minmax
fm.pmax2 <- function(o1, o2)
{
	fm.ele.wise.op(fm.bo.max, o1, o2)
}

`+.fm` <- function(o1, o2)
{
	fm.ele.wise.op(fm.bo.add, o1, o2)
}

`+.fmV` <- function(o1, o2)
{
	fm.ele.wise.op(fm.bo.add, o1, o2)
}

`-.fm` <- function(o1, o2)
{
	fm.ele.wise.op(fm.bo.sub, o1, o2)
}

`-.fmV` <- function(o1, o2)
{
	fm.ele.wise.op(fm.bo.sub, o1, o2)
}

`*.fm` <- function(o1, o2)
{
	fm.ele.wise.op(fm.bo.mul, o1, o2)
}

`*.fmV` <- function(o1, o2)
{
	fm.ele.wise.op(fm.bo.mul, o1, o2)
}

`/.fm` <- function(o1, o2)
{
	fm.ele.wise.op(fm.bo.div, o1, o2)
}

`/.fmV` <- function(o1, o2)
{
	fm.ele.wise.op(fm.bo.div, o1, o2)
}

`==.fm`  <- function(o1, o2)
{
	fm.ele.wise.op(fm.bo.eq, o1, o2)
}

`==.fmV`  <- function(o1, o2)
{
	fm.ele.wise.op(fm.bo.eq, o1, o2)
}

`!=.fm`  <- function(o1, o2)
{
	!fm.ele.wise.op(fm.bo.eq, o1, o2)
}

`!=.fmV`  <- function(o1, o2)
{
	!fm.ele.wise.op(fm.bo.eq, o1, o2)
}

`>.fm`  <- function(o1, o2)
{
	fm.ele.wise.op(fm.bo.gt, o1, o2)
}

`>.fmV`  <- function(o1, o2)
{
	fm.ele.wise.op(fm.bo.gt, o1, o2)
}

`>=.fm`  <- function(o1, o2)
{
	fm.ele.wise.op(fm.bo.ge, o1, o2)
}

`>=.fmV`  <- function(o1, o2)
{
	fm.ele.wise.op(fm.bo.ge, o1, o2)
}

`<=.fm`  <- function(o1, o2)
{
	!fm.ele.wise.op(fm.bo.gt, o1, o2)
}

`<=.fmV`  <- function(o1, o2)
{
	!fm.ele.wise.op(fm.bo.gt, o1, o2)
}

`<.fm`  <- function(o1, o2)
{
	!fm.ele.wise.op(fm.bo.ge, o1, o2)
}

`<.fmV`  <- function(o1, o2)
{
	!fm.ele.wise.op(fm.bo.ge, o1, o2)
}

`!.fm` <- function(o)
{
	fm.sapply(fm.buo.not, o)
}

`!.fmV` <- function(o)
{
	fm.sapply(fm.buo.not, o)
}

#' Aggregation on a FlashMatrixR vector/matrix.
#'
#' `fm.sum' returns the sum of all the values in the input vector/matrix.
#'
#' `fm.min' returns the minimum of all values in the input vector/matrix.
#'
#' `fm.max' returns the maximum of all values in the input vector/matrix.
#'
#' @param m The input vector/matrix
#' @return an integer if the input is an integer vector/matrix; a numeric
#' value if the input is a numeric vector/matrix.
#' @name fm.agg.impl
#' @author Da Zheng <dzheng5@@jhu.edu>
fm.sum <- function(m)
{
	fm.agg(fm.bo.add, m)
}

#' @name fm.agg.impl
fm.min <- function(m)
{
	fm.agg(fm.bo.min, m)
}

#' @name fm.agg.impl
fm.max <- function(m)
{
	fm.agg(fm.bo.max, m)
}

#' Miscellaneous Mathematical Functions
#'
#' `abs(x)' computes the absolute value of x,
#'
#' `sqrt(x)' computes the square root of x.
#'
#' `fm.ceil' returns the ceiling of all values in the input vector/matrix.
#'
#' `fm.floor' returns the floor of all values in the input vector/matrix.
#'
#' @param m a numeric vector or matrix.
#' @return `abs' returns a vector/matrix of the same type as the input
#' vector/matrix; `sqrt' returns a numeric vector/matrix.
#' @name misc.math
#' @author Da Zheng <dzheng5@@jhu.edu>
fm.abs <- function(m)
{
	fm.sapply(fm.buo.abs, m)
}

#' @name misc.math
fm.sqrt <- function(m)
{
	fm.sapply(fm.buo.sqrt, m)
}

#' @name misc.math
fm.ceil <- function(m)
{
	fm.sapply(fm.buo.ceil, m);
}

#' @name misc.math
fm.floor <- function(m)
{
	fm.sapply(fm.buo.floor, m);
}
