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

# This file contains the FlashR implementation of the functions in R base package.

#' Arithmetic Operators
#'
#' Perform unary or binary arithmetic operations on FlashR objects.
#'
#' @param e1,e2 One of the operands need to be a FlashR object. If one operand
#' is a matrix and the other is a vector, we perform the arithmetic operation
#' on the vector and every column of the matrix. If one operand is a scalar,
#' we perform the operation on the scalar with every element in the matrix or
#' the vector.
#'
#' There are a few exceptions. For example, the left operand of "^" has to be
#' a FlashR object and the right operand has to be a scalar.
#'
#' @name Arithmetic
#'
#' @examples
#' mat1 <- fm.runif.matrix(100, 10)
#' mat2 <- fm.runif.matrix(100, 10)
#' mat <- -mat1
#' mat <- mat1 + mat2
#' mat <- mat1 - mat2
#' mat <- mat1 * mat2
#' mat <- mat1 / mat2
#' mat <- mat1 ^ 2
NULL

#' @rdname Arithmetic
setMethod("-", signature(e1 = "fm", e2 = "missing"), function(e1, e2)
		  .sapply.fm(e1, fm.buo.neg))
#' @rdname Arithmetic
setMethod("-", signature(e1 = "fmV", e2 = "missing"), function(e1, e2)
		  .sapply.fmV(e1, fm.buo.neg))

#' @rdname Arithmetic
setMethod("+", signature(e1 = "fm", e2 = "fm"), function(e1, e2)
		  .mapply2.fm(e1, e2, fm.bo.add))
#' @rdname Arithmetic
setMethod("+", signature(e1 = "fmV", e2 = "fmV"), function(e1, e2)
		  .mapply2.fmV(e1, e2, fm.bo.add))

#' @rdname Arithmetic
setMethod("-", signature(e1 = "fm", e2 = "fm"), function(e1, e2)
		  .mapply2.fm(e1, e2, fm.bo.sub))
#' @rdname Arithmetic
setMethod("-", signature(e1 = "fmV", e2 = "fmV"), function(e1, e2)
		  .mapply2.fmV(e1, e2, fm.bo.sub))
#' @rdname Arithmetic
setMethod("-", signature(e1 = "fmV", e2 = "fm"), function(e1, e2)
		  -fm.mapply.col(e2, e1, fm.bo.sub))

#' @rdname Arithmetic
setMethod("*", signature(e1 = "fm", e2 = "fm"), function(e1, e2)
		  .mapply2.fm(e1, e2, fm.bo.mul))
#' @rdname Arithmetic
setMethod("*", signature(e1 = "fmV", e2 = "fmV"), function(e1, e2)
		  .mapply2.fmV(e1, e2, fm.bo.mul))

#' @rdname Arithmetic
setMethod("/", signature(e1 = "fm", e2 = "fm"), function(e1, e2)
		  .mapply2.fm(e1, e2, fm.bo.div))
#' @rdname Arithmetic
setMethod("/", signature(e1 = "fmV", e2 = "fmV"), function(e1, e2)
		  .mapply2.fmV(e1, e2, fm.bo.div))
#' @rdname Arithmetic
setMethod("/", signature(e1 = "fmV", e2 = "fm"), function(e1, e2) {
			  e1 <- fm.matrix(e1, nrow(e2), ncol(e2))
			  .mapply2.fm(e1, e2, fm.bo.div)
		  })

#' @rdname Arithmetic
setMethod("%%", signature(e1 = "fm", e2 = "fm"), function(e1, e2)
		  .mapply2.fm(e1, e2, fm.bo.mod))
#' @rdname Arithmetic
setMethod("%%", signature(e1 = "fmV", e2 = "fmV"), function(e1, e2)
		  .mapply2.fmV(e1, e2, fm.bo.mod))
#' @rdname Arithmetic
setMethod("%%", signature(e1 = "fmV", e2 = "fm"), function(e1, e2) {
			  e1 <- fm.matrix(e1, nrow(e2), ncol(e2))
			  .mapply2.fm(e1, e2, fm.bo.mod)
		  })

#' @rdname Arithmetic
setMethod("%/%", signature(e1 = "fm", e2 = "fm"), function(e1, e2)
		  .mapply2.fm(e1, e2, fm.bo.idiv))
#' @rdname Arithmetic
setMethod("%/%", signature(e1 = "fmV", e2 = "fmV"), function(e1, e2)
		  .mapply2.fmV(e1, e2, fm.bo.idiv))
#' @rdname Arithmetic
setMethod("%/%", signature(e1 = "fmV", e2 = "fm"), function(e1, e2) {
			  e1 <- fm.matrix(e1, nrow(e2), ncol(e2))
			  .mapply2.fm(e1, e2, fm.bo.idiv)
		  })

#' @rdname Arithmetic
setMethod("^", signature(e1 = "fm", e2 = "fm"), function(e1, e2)
		  .mapply2.fm(as.numeric(e1), as.numeric(e2), fm.bo.pow))
#' @rdname Arithmetic
setMethod("^", signature(e1 = "fmV", e2 = "fmV"), function(e1, e2)
		  .mapply2.fmV(as.numeric(e1), as.numeric(e2), fm.bo.pow))
#' @rdname Arithmetic
setMethod("^", signature(e1 = "fmV", e2 = "fm"), function(e1, e2) {
		  e1 <- fm.matrix(e1, nrow(e2), ncol(e2))
		  .mapply2.fm(as.numeric(e1), as.numeric(e2), fm.bo.pow)
		  })

#' Matrix multiplication
#'
#' Multiplies two matrices, if they are conformable. If one argument is
#' a vector, it will be promoted to either a row or column matrix to make
#' the two arguments conformable.
#'
#' @param x a FlashR matrix.
#' @param y can be a FlashR vector or matrix, an R vector or matrix.
#' @return a FlashR matrix.
#' @name matmult
#'
#' @examples
#' mat1 <- fm.runif.matrix(100, 10)
#' mat2 <- fm.runif.matrix(10, 10)
#' mat <- mat1 %*% mat2
NULL

#' @rdname matmult
setMethod("%*%", signature(x = "fm", y = "fm"), function(x, y) fm.multiply(x, y))
#' @rdname matmult
setMethod("%*%", signature(x = "fm", y = "fmV"), function(x, y)
		  fm.as.matrix(fm.multiply(x, y)))
#' @rdname matmult
setMethod("%*%", signature(x = "fmV", y = "fm"), function(x, y) {
		  res <- fm.as.matrix(fm.multiply(t(y), x))
		  t(res)
})
#' @rdname matmult
setMethod("%*%", signature(x = "fm", y = "ANY"),
		  function(x, y) fm.as.matrix(fm.multiply(x, fm.conv.R2FM(y))))
#' @rdname matmult
setMethod("%*%", signature(x = "ANY", y = "fm"), function(x, y) {
		  if (is.vector(x))
			  t(fm.multiply(t(y), fm.conv.R2FM(x)))
		  else
			  t(fm.multiply(t(y), t(fm.conv.R2FM(x))))
		  })

#' Relational Operators
#'
#' Binary operators which allow the comparison of values in atomic
#' vectors.
#'
#' @param e1,e2 One of the operands need to be a FlashR object. If one operand
#' is a matrix and the other is a vector, we perform the arithmetic operation
#' on the vector and every column of the matrix. If one operand is a scalar,
#' we perform the operation on the scalar with every element in the matrix or
#' the vector.
#'
#' @name Comparison
#'
#' @examples
#' mat1 <- fm.runif.matrix(100, 10)
#' mat2 <- fm.runif.matrix(100, 10)
#' mat <- mat1 == mat2
#' mat <- mat1 != mat2
#' mat <- mat1 < mat2
#' mat <- mat1 <= mat2
#' mat <- mat1 > mat2
#' mat <- mat1 >= mat2
NULL

#' @rdname Comparison
setMethod("==", signature(e1 = "fm", e2 = "fm"), function(e1, e2)
		  .mapply2.fm(e1, e2, fm.bo.eq))
#' @rdname Comparison
setMethod("==", signature(e1 = "fmV", e2 = "fmV"), function(e1, e2)
		  .mapply2.fmV(e1, e2, fm.bo.eq))

#' @rdname Comparison
setMethod("!=", signature(e1 = "fm", e2 = "fm"), function(e1, e2)
		  .mapply2.fm(e1, e2, fm.bo.neq))
#' @rdname Comparison
setMethod("!=", signature(e1 = "fmV", e2 = "fmV"), function(e1, e2)
		  .mapply2.fmV(e1, e2, fm.bo.neq))

#' @rdname Comparison
setMethod(">", signature(e1 = "fm", e2 = "fm"), function(e1, e2)
		  .mapply2.fm(e1, e2, fm.bo.gt))
#' @rdname Comparison
setMethod(">", signature(e1 = "fmV", e2 = "fmV"), function(e1, e2)
		  .mapply2.fmV(e1, e2, fm.bo.gt))
#' @rdname Comparison
setMethod(">", signature(e1 = "fmV", e2 = "fm"), function(e1, e2)
		  fm.mapply.col(e2, e1, fm.bo.lt))

#' @rdname Comparison
setMethod(">=", signature(e1 = "fm", e2 = "fm"), function(e1, e2)
		  .mapply2.fm(e1, e2, fm.bo.ge))
#' @rdname Comparison
setMethod(">=", signature(e1 = "fmV", e2 = "fmV"), function(e1, e2)
		  .mapply2.fmV(e1, e2, fm.bo.ge))
#' @rdname Comparison
setMethod(">=", signature(e1 = "fmV", e2 = "fm"), function(e1, e2)
		  fm.mapply.col(e2, e1, fm.bo.le))

#' @rdname Comparison
setMethod("<=", signature(e1 = "fm", e2 = "fm"), function(e1, e2)
		  .mapply2.fm(e1, e2, fm.bo.le))
#' @rdname Comparison
setMethod("<=", signature(e1 = "fmV", e2 = "fmV"), function(e1, e2)
		  .mapply2.fmV(e1, e2, fm.bo.le))
#' @rdname Comparison
setMethod("<=", signature(e1 = "fmV", e2 = "fm"), function(e1, e2)
		  fm.mapply.col(e2, e1, fm.bo.ge))

#' @rdname Comparison
setMethod("<", signature(e1 = "fm", e2 = "fm"), function(e1, e2)
		  .mapply2.fm(e1, e2, fm.bo.lt))
#' @rdname Comparison
setMethod("<", signature(e1 = "fmV", e2 = "fmV"), function(e1, e2)
		  .mapply2.fmV(e1, e2, fm.bo.lt))
#' @rdname Comparison
setMethod("<", signature(e1 = "fmV", e2 = "fm"), function(e1, e2)
		  fm.mapply.col(e2, e1, fm.bo.gt))

#' Logical Operators
#'
#' These operators act on logical and number-like vectors.
#'
#' \code{!} indicates logical negation (NOT).
#' \code{&} indicates logical AND and \code{|} indicates logical OR.
#'
#' @param e1,e2 One of the operands need to be a FlashR object. If one operand
#' is a matrix and the other is a vector, we perform the arithmetic operation
#' on the vector and every column of the matrix. If one operand is a scalar,
#' we perform the operation on the scalar with every element in the matrix or
#' the vector.
#'
#' @name Logic
#'
#' @examples
#' tmp1 <- fm.runif.matrix(100, 10)
#' tmp2 <- fm.runif.matrix(100, 10)
#' mat1 <- tmp1 > tmp2
#' mat2 <- tmp1 < tmp2
#' mat1 | mat2
#' mat1 & mat2
#' !mat1
NULL

#' @rdname Logic
setMethod("|", signature(e1 = "fm", e2 = "fm"), function(e1, e2)
		  .mapply2.fm(as.logical(e1), as.logical(e2), fm.bo.or))
#' @rdname Logic
setMethod("|", signature(e1 = "fmV", e2 = "fmV"), function(e1, e2)
		  .mapply2.fmV(as.logical(e1), as.logical(e2), fm.bo.or))

#' @rdname Logic
setMethod("&", signature(e1 = "fm", e2 = "fm"), function(e1, e2)
		  .mapply2.fm(as.logical(e1), as.logical(e2), fm.bo.and))
#' @rdname Logic
setMethod("&", signature(e1 = "fmV", e2 = "fmV"), function(e1, e2)
		  .mapply2.fmV(as.logical(e1), as.logical(e2), fm.bo.and))

#' @rdname Logic
`!.fm` <- function(e1)
{
	.sapply.fm(as.logical(e1), fm.buo.not)
}

#' @rdname Logic
`!.fmV` <- function(e1)
{
	.sapply.fmV(as.logical(e1), fm.buo.not)
}

#' Maxima and Minima
#'
#' Returns the (parallel) maxima and minima of the input values.
#'
#' pmax2 returns the maximum of two objects.
#' pmin2 returns the minimum of two objects.
#'
#' If one of the arguments is a matrix and the other is a vector, \code{pmax2}
#' or \code{pmin2} return a matrix, each of whose column is parallel maxima
#' and minima of the column of the input matrix and the input vector.
#'
#' @param e1,e2 One of the operands need to be a FlashR object. If one
#' operand is a matrix and the other is a vector, we perform the operation
#' on the vector and every column of the matrix. If one operand is a scalar,
#' we perform the operation on the scalar with every element in the matrix or
#' the vector.
#' @param x a FlashR vector or matrix.
#' @param ... FlashR vectors or matrices.
#' @param na.rm a logical indicating whether missing values should be removed.
#' @return \code{pmin2}, \code{pmax2}, \code{pmin} and \code{pmax} return
#'         a FlashR object, and \code{min} and \code{max} return an R scalar.
#'
#' @name Extremes
#'
#' @examples
#' mat1 <- fm.runif.matrix(100, 10)
#' mat2 <- fm.runif.matrix(100, 10)
#' res <- pmax2(mat1, mat2)
#' res <- pmin2(mat1, mat2)
#' res <- pmax(mat1, mat2)
#' res <- pmin(mat1, mat2)
#' res <- max(mat1, mat2)
#' res <- min(mat1, mat2)
pmax2 <- function(e1, e2) {
	# if e1 is a vector and e2 is a matrix, pmax(e1, e2) returns a vector,
	# which is inconsistent with others, so in this case, we need to switch
	# the two arguments.
	if (is.vector(e1) && is.matrix(e2))
		pmax(e2, e1)
	else
		pmax(e1, e2)
}
#' @rdname Extremes
pmin2 <- function(e1, e2) {
	if (is.vector(e1) && is.matrix(e2))
		pmin(e2, e1)
	else
		pmin(e1, e2)
}
setGeneric("pmax2")
setGeneric("pmin2")

#' @rdname Extremes
setMethod("pmax2", signature(e1 = "fm", e2 = "fm"), function(e1, e2)
		  .mapply2.fm(e1, e2, fm.bo.max))
#' @rdname Extremes
setMethod("pmax2", signature(e1 = "fmV", e2 = "fmV"), function(e1, e2)
		  .mapply2.fmV(e1, e2, fm.bo.max))
#' @rdname Extremes
setMethod("pmax2", signature(e1 = "fm", e2 = "fmV"), function(e1, e2)
		  fm.mapply.col(e1, e2, fm.bo.max))
#' @rdname Extremes
setMethod("pmax2", signature(e1 = "fmV", e2 = "fm"), function(e1, e2)
		  fm.mapply.col(e2, e1, fm.bo.max))
#' @rdname Extremes
setMethod("pmax2", signature(e1 = "fm", e2 = "matrix"), function(e1, e2)
		  .mapply2.fm(e1, fm.as.matrix(e2), fm.bo.max))
#' @rdname Extremes
setMethod("pmax2", signature(e1 = "matrix", e2 = "fm"), function(e1, e2)
		  .mapply2.fm(fm.as.matrix(e1), e2, fm.bo.max))
#' @rdname Extremes
setMethod("pmax2", signature(e1 = "fm", e2 = "ANY"), function(e1, e2)
		  .mapply2.fm(e1, fm.as.matrix(e2), fm.bo.max))
#' @rdname Extremes
setMethod("pmax2", signature(e1 = "ANY", e2 = "fm"), function(e1, e2)
		  .mapply2.fm(e2, fm.as.matrix(e1), fm.bo.max))
#' @rdname Extremes
setMethod("pmax2", signature(e1 = "fmV", e2 = "ANY"), function(e1, e2)
		  pmax2(e1, fm.conv.R2FM(e2)))
#' @rdname Extremes
setMethod("pmax2", signature(e1 = "ANY", e2 = "fmV"), function(e1, e2)
		  pmax2(fm.conv.R2FM(e1), e2))

#' @rdname Extremes
setMethod("pmin2", signature(e1 = "fm", e2 = "fm"), function(e1, e2)
		  .mapply2.fm(e1, e2, fm.bo.min))
#' @rdname Extremes
setMethod("pmin2", signature(e1 = "fmV", e2 = "fmV"), function(e1, e2)
		  .mapply2.fmV(e1, e2, fm.bo.min))
#' @rdname Extremes
setMethod("pmin2", signature(e1 = "fm", e2 = "fmV"), function(e1, e2)
		  fm.mapply.col(e1, e2, fm.bo.min))
#' @rdname Extremes
setMethod("pmin2", signature(e1 = "fmV", e2 = "fm"), function(e1, e2)
		  fm.mapply.col(e2, e1, fm.bo.min))
#' @rdname Extremes
setMethod("pmin2", signature(e1 = "fm", e2 = "matrix"), function(e1, e2)
		  .mapply2.fm(e1, fm.as.matrix(e2), fm.bo.min))
#' @rdname Extremes
setMethod("pmin2", signature(e1 = "matrix", e2 = "fm"), function(e1, e2)
		  .mapply2.fm(fm.as.matrix(e1), e2, fm.bo.min))
#' @rdname Extremes
setMethod("pmin2", signature(e1 = "fm", e2 = "ANY"), function(e1, e2)
		  .mapply2.fm(e1, fm.as.matrix(e2), fm.bo.min))
#' @rdname Extremes
setMethod("pmin2", signature(e1 = "ANY", e2 = "fm"), function(e1, e2)
		  .mapply2.fm(e2, fm.as.matrix(e1), fm.bo.min))
#' @rdname Extremes
setMethod("pmin2", signature(e1 = "fmV", e2 = "ANY"), function(e1, e2)
		  pmin2(e1, fm.conv.R2FM(e2)))
#' @rdname Extremes
setMethod("pmin2", signature(e1 = "ANY", e2 = "fmV"), function(e1, e2)
		  pmin2(fm.conv.R2FM(e1), e2))

#' @rdname Arithmetic
setMethod("Ops", signature(e1 = "fm", e2 = "fmV"), function(e1, e2)
		  callGeneric(e1, fm.as.matrix(e2)))
#' @rdname Arithmetic
setMethod("Ops", signature(e1 = "fmV", e2 = "fm"), function(e1, e2)
		  callGeneric(fm.as.matrix(e1), e2))
#' @rdname Arithmetic
setMethod("Ops", signature(e1 = "fm", e2 = "matrix"), function(e1, e2)
		  callGeneric(e1, fm.as.matrix(e2)))
#' @rdname Arithmetic
setMethod("Ops", signature(e1 = "matrix", e2 = "fm"), function(e1, e2)
		  callGeneric(fm.as.matrix(e1), e2))
#' @rdname Arithmetic
setMethod("Ops", signature(e1 = "fm", e2 = "ANY"), function(e1, e2)
		  callGeneric(e1, fm.as.matrix(e2)))
#' @rdname Arithmetic
setMethod("Ops", signature(e1 = "ANY", e2 = "fm"), function(e1, e2)
		  callGeneric(fm.as.matrix(e1), e2))
#' @rdname Arithmetic
setMethod("Ops", signature(e1 = "fmV", e2 = "ANY"), function(e1, e2)
		  callGeneric(e1, fm.conv.R2FM(e2)))
#' @rdname Arithmetic
setMethod("Ops", signature(e1 = "ANY", e2 = "fmV"), function(e1, e2)
		  callGeneric(fm.conv.R2FM(e1), e2))

.min.int <- function(x, ..., na.rm) {
	others <- list(...)
	if (na.rm) {
		max.val <- .get.max.val(typeof(x))
		x <- .replace.na(x, max.val)
		others <- .replace.na.list(others, max.val)
	}
	res <- fm.agg(x, fm.bo.min)
	if (length(others) >= 1) {
		res <- list(res)
		for (arg in others)
			res <- c(res, fm.agg(arg, fm.bo.min))
		res <- lapply(res, function(o) fm.conv.FM2R(o))
		res <- min(unlist(res))
	}
	res
}

.max.int <- function(x, ..., na.rm) {
	others <- list(...)
	if (na.rm) {
		min.val <- .get.min.val(typeof(x))
		x <- .replace.na(x, min.val)
		others <- .replace.na.list(others, min.val)
	}
	res <- fm.agg(x, fm.bo.max)
	if (length(others) >= 1) {
		res <- list(res)
		for (arg in others)
			res <- c(res, fm.agg(arg, fm.bo.max))
		res <- lapply(res, function(o) fm.conv.FM2R(o))
		res <- max(unlist(res))
	}
	else
		res
}

.pmin.int <- function(..., na.rm = FALSE)
{
	args <- list(...)
	if (length(args) == 0)
		stop("no arguments")
	else if (length(args) == 1)
		return(args[[1]])

	if (na.rm) {
		args <- .replace.na.list(args, .get.max.val(typeof(args[[1]])))
		.mapply.list(args, fm.bo.min)
		# TODO there is a bug for the case that all elements in a location
		# is NA.
	}
	else
		.mapply.list(args, fm.bo.min)
}

.pmax.int <- function(..., na.rm = FALSE)
{
	args <- list(...)
	if (length(args) == 0)
		stop("no arguments")
	else if (length(args) == 1)
		return(args[[1]])

	if (na.rm) {
		args <- .replace.na.list(args, .get.min.val(typeof(args[[1]])))
		.mapply.list(args, fm.bo.max)
		# TODO there is a bug for the case that all elements in a location
		# is NA.
	}
	else
		.mapply.list(args, fm.bo.max)
}

#' @rdname Extremes
setMethod("min", "fm", .min.int)
#' @rdname Extremes
setMethod("min", "fmV", .min.int)

#' @rdname Extremes
setMethod("max", "fm", .max.int)
#' @rdname Extremes
setMethod("max", "fmV", .max.int)

setGeneric("pmax", signature="...")
setGeneric("pmin", signature="...")

#' @rdname Extremes
setMethod("pmin", "fm", .pmin.int)
#' @rdname Extremes
setMethod("pmin", "fmV", .pmin.int)

#' @rdname Extremes
setMethod("pmax", "fm", .pmax.int)
#' @rdname Extremes
setMethod("pmax", "fmV", .pmax.int)

#' Dimensions of an FlashR Object
#'
#' Retrieve the dimension of an FlashR object.
#'
#' @param x A FlashR matrix.
#'
#' @name dim
#'
#' @examples
#' mat <- fm.runif.matrix(100, 10)
#' dim(mat)
NULL

#' @rdname dim
setMethod("dim", signature(x = "fm"), function(x) c(x@nrow, x@ncol))

#' The Number of Rows/Columns of a matrix.
#'
#' \code{nrow} and \code{ncol} return the number of rows or columns
#' present in \code{x}.
#'
#' @param x a matrix.
#' @return an integer of length 1 or NULL.
#' @name nrow
#'
#' @examples
#' mat <- fm.runif.matrix(100, 10)
#' nrow(mat)
#' ncol(mat)
NULL

#' @rdname nrow
setMethod("nrow", signature(x = "fm"), function(x) x@nrow)
#' @rdname nrow
setMethod("ncol", signature(x = "fm"), function(x) x@ncol)

#' Length of an Object
#'
#' Get the length of a vector or a matrix.
#'
#' @param x an FlashR object.
#' @return an integer of length 1 or a double if the object has more than
#' 2^31-1 elements.
#' @name length
#'
#' @examples
#' mat <- fm.runif.matrix(100, 10)
#' length(mat)
NULL

#' @rdname length
setMethod("length", signature(x = "fmV"), function(x) x@len)
#' @rdname length
setMethod("length", signature(x = "fm"), function(x) x@nrow * x@ncol)

#' The Type of an Object
#'
#' \code{typeof} determines the R type of an FlashR object.
#'
#' @param x a FlashR object.
#' @return A character string. Current values are "logical", "integer",
#' "double".
#' @name typeof
#'
#' @examples
#' mat <- fm.runif.matrix(100, 10)
#' typeof(mat)
NULL

#' @rdname typeof
setMethod("typeof", signature(x = "fm"), function(x) .typeof.int(x))
#' @rdname typeof
setMethod("typeof", signature(x = "fmV"), function(x) .typeof.int(x))

#' Dimnames of an Object
#'
#' Retrieve or set the dimnames of a matrix.
#'
#' @param x a FlashR matrix.
#' @param value a list that provides values for dimension names.
#' @return \code{dimnames} return NULL if there aren't dimension names.
#'         Otherwise, return the dimension names.
#' @name dimnames
#'
#' @examples
#' mat <- fm.runif.matrix(100, 10)
#' dimnames(mat)
#' dimnames(mat) <- list("dim1", "dim2")
NULL

.get.dimnames <- function(x)
{
	if (!is.null(x@attrs[["dimnames"]])) x@attrs$dimnames else NULL
}

.set.dimnames <- function(x, value)
{
	if (is.null(x@attrs))
		x@attrs <- list()
	if (is.null(x@attrs[["dimnames"]]))
		x@attrs$dimnames <- list(NULL, NULL)
	if (length(value) != length(dim(x)))
		stop("length of dimnames not equal to array extent")
	for (i in 1:length(value))
		x@attrs$dimnames[[i]] <- as.character(value[[i]])
	x
}

#' @rdname dimnames
setMethod("dimnames", signature(x = "fm"), .get.dimnames)
#' @rdname dimnames
setMethod("dimnames", signature(x = "fmV"), .get.dimnames)

#' @rdname dimnames
setMethod("dimnames<-", signature(x = "fm", value="list"), .set.dimnames)
#' @rdname dimnames
setMethod("dimnames<-", signature(x = "fmV", value="list"), .set.dimnames)

#' The Names of an Object
#'
#' Functions to get or set the names of an object.
#'
#' Currently, \code{value} has to have the same length as \code{x} and has
#' the same shape. FlashR currently doesn't support character vectors yet.
#'
#' @param x a FlashR array
#' @param value a FlashR array.
#' @return \code{names} returns \code{NULL} or a FlashR array of the same
#' length as \code{x}. \code{names<-} returns the updated object.
#' @name names
NULL

.get.names <- function(x)
{
	if (!is.null(x@attrs[["names"]])) x@attrs$names else NULL
}

.set.names <- function(x, value)
{
	if (length(x) != length(value))
		stop("The number of names isn't equal to the number of elements")
	# if two inputs have different dimensions.
	if (length(dim(x)) != length(dim(value))
		# if both are matrices but they have different shapes.
		|| (!is.null(dim(x)) && any(dim(x) != dim(value))))
		stop("x and value have different shapes.")
	if (is.null(x@attrs))
		x@attrs <- list()
	# TODO I should make sure the names are characters.
	x@attrs$names <- value
	x
}

#' @rdname names
setMethod("names", signature(x = "fm"), .get.names)
#' @rdname names
setMethod("names", signature(x = "fmV"), .get.names)

#' @rdname names
setMethod("names<-", signature(x = "fm", value="fm"), .set.names)
#' @rdname names
setMethod("names<-", signature(x = "fmV", value="fmV"), .set.names)

#' Matrix Transpose
#'
#' Given a matrix \code{x}, \code{t} returns the transpose of \code{x}.
#'
#' @param x a FlashR matrix.
#' @return a matrix.
#' @name transpose
#'
#' @examples
#' mat <- fm.runif.matrix(100, 10)
#' t(mat)
NULL

#' @rdname transpose
setMethod("t", signature(x = "fm"), function(x) fm.t(x))
#' @rdname transpose
setMethod("t", signature(x = "fmV"), function(x) fm.t(fm.as.matrix(x)))

.mapply.list <- function(data, FUN)
{
	n <- length(data)
	res <- data[[1]]
	if (n > 1) {
		for (arg in data[2:n]) {
			if (class(res) == class(arg))
				res <- fm.mapply2(res, arg, FUN)
			else if (class(res) == "fm" && class(arg) == "fmV")
				res <- fm.mapply.col(res, arg, FUN)
			# We assume that FUN is commutative.
			# This function is used for pmin/pmax, so it's fine.
			else if (class(res) == "fmV" && class(arg) == "fm")
				res <- fm.mapply.col(arg, res, FUN)
			else
				stop("unknown arguments")
		}
	}
	res
}

.get.zero <- function(type)
{
	if (type == "double")
		0
	else if (type == "integer")
		as.integer(0)
	else if (type == "logical")
		FALSE
	else
		stop("only integer and numeric value have 0")
}

.get.max.val <- function(type)
{
	if (type == "double")
		.Machine$double.xmax
	else if (type == "integer")
		.Machine$integer.max
	else
		TRUE
}

.get.min.val <- function(type)
{
	if (type == "double")
		-.Machine$double.xmax
	else if (type == "integer")
		-.Machine$integer.max
	else
		FALSE
}

.range1 <- function(x, na.rm)
{
	if (typeof(x) == "logical")
		x <- as.integer(x)
	if (na.rm) {
		x.min <- ifelse(is.na(x), .get.min.val(typeof(x)), x)
		x.max <- ifelse(is.na(x), .get.max.val(typeof(x)), x)
		tmp1 <- fm.agg(x.max, fm.bo.min)
		tmp2 <- fm.agg(x.min, fm.bo.max)
	}
	else {
		tmp1 <- fm.agg(x, fm.bo.min)
		tmp2 <- fm.agg(x, fm.bo.max)
	}
	c(.fmV2scalar(tmp1), .fmV2scalar(tmp2))
}

# We replace both NA and NaN
.replace.na <- function(obj, val)
{
	ifelse(is.na(obj), val, obj)
}
.replace.na.list <- function(objs, val)
{
	if (length(objs) == 0)
		return(objs)

	for (i in 1:length(objs))
		objs[[i]] <- ifelse(is.na(objs[[i]]), val, objs[[i]])
	objs
}

#' Are some Values True?
#'
#' \code{any} tests whether some of the values in a set of logical vectors are
#' true? \code{fm.any} tests whether some of the values in a single vector are
#' true. \code{fm.any} can evaluate it lazily.
#'
#' @param x a logical FlashR vector.
#' @param ... zero or more logical vectors.
#' @param lazy indicates whether or not to evaluate it lazily.
#' @param na.rm a logical indicating whether missing values should be removed.
#' @return The value is a logical vector of length one.
#' @name any
#'
#' @examples
#' mat <- fm.runif.matrix(100, 10) > 0.5
#' any(mat)
#' fm.any(mat)
fm.any <- function(x, lazy=FALSE)
{
	ret <- fm.agg(x, fm.bo.or)
	if (!lazy)
		.fmV2scalar(ret)
	else
		ret
}

.any.int <- function(x, ..., na.rm)
{
	others <- list(...)
	# If the inputs aren't logical, we should cast them
	# to logical values.
	if (typeof(x) != "logical")
		x <- as.logical(x)
	if (length(others) > 0) {
		for (i in 1:length(others)) {
			if (typeof(others[[i]]) != "logical")
				others[[i]] <- as.logical(others[[i]])
		}
	}
	if (na.rm) {
		x <- .replace.na(x, FALSE)
		others <- .replace.na.list(others, FALSE)
	}
	res <- fm.conv.FM2R(fm.agg(x, fm.bo.or))
	if (length(others) >= 1) {
		if (!is.na(res) && res)
			return(TRUE)
		for (arg in others) {
			tmp <- fm.conv.FM2R(fm.agg(arg, fm.bo.or))
			if (!is.na(tmp) && tmp)
				return(TRUE)
			res <- res | tmp
		}
	}
	res
}

#' @rdname any
setMethod("any", "fm", .any.int)
#' @rdname any
setMethod("any", "fmV", .any.int)

#' Are All Values True?
#'
#' \code{all} tests whether all of the values in a set of logical vectors are
#' true? \code{fm.all} tests whether all values in a single vector are true.
#' \code{fm.all} can evaluate it lazily.
#'
#' @param x a logical FlashR vector.
#' @param ... zero or more logical vectors.
#' @param lazy indicates whether or not to evaluate it lazily.
#' @param na.rm a logical indicating whether missing values should be removed.
#' @return The value is a logical vector of length one.
#' @name all
#'
#' @examples
#' mat <- fm.runif.matrix(100, 10) > 0.5
#' all(mat)
#' fm.all(mat)
fm.all <- function(x, lazy=FALSE)
{
	ret <- fm.agg(x, fm.bo.and)
	if (!lazy)
		.fmV2scalar(ret)
	else
		ret
}

.all.int <- function(x, ..., na.rm)
{
	others <- list(...)
	# If the inputs aren't logical, we should cast them
	# to logical values.
	if (typeof(x) != "logical")
		x <- as.logical(x)
	if (length(others) > 0) {
		for (i in 1:length(others)) {
			if (typeof(others[[i]]) != "logical")
				others[[i]] <- as.logical(others[[i]])
		}
	}
	if (na.rm) {
		x <- .replace.na(x, TRUE)
		others <- .replace.na.list(others, TRUE)
	}
	res <- fm.conv.FM2R(fm.agg(x, fm.bo.and))
	if (length(others) >= 1) {
		if (!is.na(res) && !res)
			return(FALSE)
		for (arg in others) {
			tmp <- fm.conv.FM2R(fm.agg(arg, fm.bo.and))
			if (!is.na(tmp) && !tmp)
				return(FALSE)
			res <- res & tmp
		}
	}
	res
}

#' @rdname all
setMethod("all", "fm", .all.int)
#' @rdname all
setMethod("all", "fmV", .all.int)

#' Sum of Vector Elements
#'
#' \code{sum} returns the sum of all the values in its arguments.
#'
#' @param x a FlashR vector or matrix.
#' @param ... zero or more vectors or matrices.
#' @param lazy indicates whether or not to evaluate it lazily.
#' @param na.rm logical. Should missing values (including \code{NaN}) be removed?
#' @return The sum.
#' @name sum
#'
#' @examples
#' mat <- fm.runif.matrix(100, 10)
#' sum(mat)
NULL

.sum.int <- function(x, ..., na.rm)
{
	# TODO we need to handle na.rm for all of the functions here
	# properly.
	others <- list(...)
	if (na.rm) {
		zero <- .get.zero(typeof(x))
		x <- .replace.na(x, zero)
		others <- .replace.na.list(others, zero)
	}
	res <- fm.agg(x, fm.bo.add)
	if (length(others) >= 1) {
		res <- list(res)
		for (arg in others)
			res <- c(res, fm.agg(arg, fm.bo.add))
		res <- lapply(res, function(o) fm.conv.FM2R(o))
		res <- sum(unlist(res))
	}
	res
}

#' @rdname sum
setMethod("sum", "fm", .sum.int)
#' @rdname sum
setMethod("sum", "fmV", .sum.int)

.range.int <- function(x, ..., na.rm)
{
	others <- list(...)
	res <- .range1(x, na.rm)
	if (length(others) >= 1) {
		for (arg in others) {
			tmp <- .range1(arg, na.rm)
			# If there is NA in the data, there is no point of
			# continuing.
			if (is.na(tmp))
				return(tmp)
			res[1] <- min(res[1], tmp[1])
			res[2] <- max(res[2], tmp[2])
		}
	}
	res
}

#' Range of Values
#'
#' \code{range} returns a vector containing the minimum and maximum of all
#' the given arguments.
#'
#' @param x a FlashR vector or matrix.
#' @param ... any FlashR vectors or matrices.
#' @param na.rm logical, indicating if \code{NA} should be omitted.
#' @name range
#'
#' @examples
#' mat <- fm.runif.matrix(100, 10)
#' range(mat)
NULL

#' @rdname range
setMethod("range", "fm", .range.int)
#' @rdname range
setMethod("range", "fmV", .range.int)

.mean.int <- function(x, ...)
{
	# TODO I need to implemented trimmed mean
	args <- list(...)
	na.rm <- FALSE
	if (!is.null(args[["na.rm"]]))
		na.rm <- args[["na.rm"]]
	n <- length(x)
	# TODO I need to lazily calculate it.
	if (na.rm)
		n <- n - sum(is.na(x))
	sum(as.double(x), na.rm=na.rm) / n
}

#' Arithmetic Mean
#'
#' Compute arithmetic mean.
#'
#' @param x A FlashR vector or matrix.
#' @param ... further arguments passed to or from other methods.
#' @name mean
#'
#' @examples
#' mat <- fm.runif.matrix(100, 10)
#' mean(mat)
NULL

#' @rdname mean
setMethod("mean", "fm", .mean.int)
#' @rdname mean
setMethod("mean", "fmV", .mean.int)

#' Miscellaneous Mathematical Functions
#'
#' \code{abs(x)} computes the absolute value of x, \code{sqrt(x)} computes the square
#' root of x.
#'
#' @param x a numeric vector or array.
#' @name MathFun
#'
#' @examples
#' mat <- abs(fm.runif.matrix(100, 10))
#' mat <- sqrt(fm.runif.matrix(100, 10))
NULL

#' @rdname MathFun
setMethod("abs", signature(x = "fm"), function(x) .sapply.fm(x, fm.buo.abs))
#' @rdname MathFun
setMethod("abs", signature(x = "fmV"), function(x) .sapply.fmV(x, fm.buo.abs))
#' @rdname MathFun
setMethod("sqrt", signature(x = "fm"), function(x) .sapply.fm(x, fm.buo.sqrt))
#' @rdname MathFun
setMethod("sqrt", signature(x = "fmV"), function(x) .sapply.fmV(x, fm.buo.sqrt))

#' Rounding of Numbers
#'
#' \code{ceiling} takes a single numeric argument \code{x} and returns
#' a numeric vector containing the smallest integers not less than
#' the corresponding elements of \code{x}.
#'
#' \code{floor} takes a single numeric argument \code{x} and returns a numeric
#' vector containing the largest integers not greater than the corresponding
#' elements of \code{x}.
#'
#' \code{round} rounds the values in its first argument to the specified number
#' of decimal places (default 0).
#'
#' @param x a numeric vector.
#' @name round
#'
#' @examples
#' mat <- ceiling(fm.runif.matrix(100, 10))
#' mat <- floor(fm.runif.matrix(100, 10))
#' mat <- round(fm.runif.matrix(100, 10))
NULL

#' @rdname round
setMethod("ceiling", signature(x = "fm"), function(x) .sapply.fm(as.numeric(x), fm.buo.ceil))
#' @rdname round
setMethod("ceiling", signature(x = "fmV"), function(x) .sapply.fmV(as.numeric(x), fm.buo.ceil))
#' @rdname round
setMethod("floor", signature(x = "fm"), function(x) .sapply.fm(as.numeric(x), fm.buo.floor))
#' @rdname round
setMethod("floor", signature(x = "fmV"), function(x) .sapply.fmV(as.numeric(x), fm.buo.floor))
#' @rdname round
setMethod("round", signature(x = "fm"), function(x) .sapply.fm(as.numeric(x), fm.buo.round))
#' @rdname round
setMethod("round", signature(x = "fmV"), function(x) .sapply.fmV(as.numeric(x), fm.buo.round))

#' Logarithms and Exponentials
#'
#' \code{log} computes logarithms, by default natural logarithms, \code{log10}
#' computes common (i.e., base 10) logarithms, and \code{log2} computes binary
#' (i.e., base 2) logarithms. The general form log(x, base) computes logarithms
#' with \code{base}.
#'
#' \code{exp} computes the exponential function.
#'
#' @param x a numeric vector.
#' @param base a positive number.
#' @return A vector of the same length as \code{x} containing the transformed
#' values. \code{log(0)} gives \code{-Inf}, and \code{log(x)} for negative
#' values of \code{x} is \code{NaN}.  \code{exp(-Inf)} is 0.
#' @name log
#'
#' @examples
#' mat <- log(fm.runif.matrix(100, 10))
#' mat <- log10(fm.runif.matrix(100, 10))
#' mat <- log2(fm.runif.matrix(100, 10))
#' mat <- exp(fm.runif.matrix(100, 10))
NULL

#' @rdname log
setMethod("log10", signature(x = "fm"), function(x) .sapply.fm(as.numeric(x), fm.buo.log10))
#' @rdname log
setMethod("log10", signature(x = "fmV"), function(x) .sapply.fmV(as.numeric(x), fm.buo.log10))
#' @rdname log
setMethod("log2", signature(x = "fm"), function(x) .sapply.fm(as.numeric(x), fm.buo.log2))
#' @rdname log
setMethod("log2", signature(x = "fmV"), function(x) .sapply.fmV(as.numeric(x), fm.buo.log2))
#' @rdname log
setMethod("exp", signature(x = "fm"), function(x) fm.mapply2(exp(1), as.numeric(x), fm.bo.pow))
#' @rdname log
setMethod("exp", signature(x = "fmV"), function(x) fm.mapply2(exp(1), as.numeric(x), fm.bo.pow))
#' @rdname log
setMethod("log", "fm", function(x, base=exp(1)) {
		  if (base == exp(1))
			  .sapply.fm(as.numeric(x), fm.buo.log)
		  else
			  .sapply.fm(as.numeric(x), fm.buo.log) / log(base)
})
#' @rdname log
setMethod("log", "fmV", function(x, base=exp(1)) {
		  if (base == exp(1))
			  .sapply.fmV(as.numeric(x), fm.buo.log)
		  else
			  .sapply.fmV(as.numeric(x), fm.buo.log) / log(base)
})

#' Form Row and Column Sums and Means
#'
#' Form row and column sums and means for numeric arrays.
#'
#' @param x a FlashR matrix.
#' @param na.rm logical. Should missing values (including NaN) be omitted
#' from the calculations?
#' @return a FlashR vector.
#' @name colSums
#'
#' @examples
#' mat <- fm.runif.matrix(100, 10)
#' mat <- rowSums(mat)
#' mat <- colSums(mat)
#' mat <- rowMeans(mat)
#' mat <- colMeans(mat)
NULL

#' @rdname colSums
setMethod("rowSums", signature(x = "fm", na.rm = "ANY"),
		  function(x, na.rm) {
			  if (na.rm)
				  x <- .replace.na(x, .get.zero(typeof(x)))
			  fm.agg.mat(x, 1, fm.bo.add)
		  })
#' @rdname colSums
setMethod("colSums", signature(x = "fm", na.rm = "ANY"),
		  function(x, na.rm) {
			  if (na.rm)
				  x <- .replace.na(x, .get.zero(typeof(x)))
			  fm.agg.mat(x, 2, fm.bo.add)
		  })
#' @rdname colSums
setMethod("rowMeans", signature(x = "fm", na.rm = "ANY"),
		  function(x, na.rm) {
			  x <- as.double(x)
			  if (na.rm) {
				  sums <- rowSums(x, na.rm)
				  nums <- rowSums(!is.na(x), FALSE)
				  sums / nums
			  }
			  else
				  rowSums(x, na.rm) / dim(x)[2]
		  })
#' @rdname colSums
setMethod("colMeans", signature(x = "fm", na.rm = "ANY"),
		  function(x, na.rm) {
			  x <- as.double(x)
			  if (na.rm) {
				  sums <- colSums(x, na.rm)
				  nums <- colSums(!is.na(x), FALSE)
				  sums / nums
			  }
			  else
				  colSums(x, na.rm) / dim(x)[1]
		  })

#' Print FlashR objects.
#'
#' \code{print} prints the information of a FlashR object.
#'
#' @param x a FlashR object. It can be a vector, a matrix or a FlashR
#' binary operator.
#' @name print
NULL

#' @rdname print
setMethod("print", signature(x = "fm"), function(x)
		  cat("FlashR matrix ", x@name, ": ", dim(x)[1], " rows, ", dim(x)[2],
			  " columns, is sparse: ", fm.is.sparse(x), "\n", sep=""))
#' @rdname print
setMethod("print", signature(x = "fmV"), function(x)
	cat("FlashR vector ", x@name, ": length: ", length(x), "\n", sep=""))
#' @rdname print
setMethod("print", signature(x = "fmVFactor"), function(x)
	cat("FlashR factor vector ", x@name, ": length: ", length(x),
		", max level: ", x@num.levels, "\n", sep=""))
#' @rdname print
setMethod("print", signature(x = "fm.bo"), function(x)
	cat("FLashR basic operator:", x@name, "\n"))

#' Count the number of elements
#'
#' \code{fm.table} counts the number of occurences of each unique value
#' in a FlashR vector.
#'
#' @param x a FlashR vector or a "table" object for \code{as.vector} and
#'          \code{as.data.frame}.
#' @return a table object.
#' \itemize{
#' \item{val}{The unique values in the vector.}
#' \item{Freq}{The number of occurences of each unique value.}
#' }
#' @name fm.table
#'
#' @examples
#' vec <- as.integer(fm.runif(1000, min=0, max=10))
#' tbl <- fm.table(vec)
#' cnts <- as.vector(tbl)
#' df <- as.data.frame(tbl)
NULL

setClass("fm.table", representation(val = "fmV", Freq = "fmV"))

#' @rdname fm.table
fm.table <- function(x)
{
	valid.vec <- function(x) x@len > 0
	if (class(x) == "fmVFactor" && valid.vec(x@vals) && valid.vec(x@cnts))
		ret <- list(val=x@vals, agg=x@cnts)
	else {
		count <- fm.create.agg.op(fm.bo.count, fm.bo.add, "count")
		ret <- fm.sgroupby(x, count)
	}
	if (!is.null(ret))
		new("fm.table", val=ret$val, Freq=ret$agg)
	else
		NULL
}

#' @rdname fm.table
setMethod("as.vector", signature(x = "fm.table"), function(x) x@Freq)
#' @rdname fm.table
setMethod("as.data.frame", signature(x = "fm.table"),
		  function(x, row.names = NULL, optional = FALSE, ...)
			  list(val=x@val, Freq=x@Freq))

#' Logical Vectors
#'
#' \code{as.logical} coerces objects of type \code{"logical"}.
#' \code{is.logical} is a more general test of an object being
#' interpretable as logicals.
#'
#' @param x a FlashR object to be coerced or tested.
#' @return \code{as.logical} returns a logical FlashR object,
#' \code{is.logical} returns a logical value.
#' @name logical
NULL

#' @rdname logical
setMethod("as.logical", "fm", function(x) {
		  if (.typeof.int(x) == "logical")
			  x
		  # Even though the underlying matrix uses the same C type to store
		  # logical values and integers, we still need to cast it to cast
		  # NA properly.
		  else
			  fm.sapply(x, fm.buo.as.logical)
	})
#' @rdname logical
setMethod("as.logical", "fmV", function(x) {
		  if (.typeof.int(x) == "logical")
			  x
		  # Even though the underlying matrix uses the same C type to store
		  # logical values and integers, we still need to cast it to cast
		  # NA properly.
		  else
			  fm.sapply(x, fm.buo.as.logical)
	})

#' Integer Vectors
#'
#' \code{as.integer} coerces objects of type \code{"integer"}.
#' \code{is.integer} is a more general test of an object being
#' interpretable as integers.
#'
#' @param x a FlashR object to be coerced or tested.
#' @return \code{as.integer} returns an integer FlashR object,
#' \code{is.integer} returns a logical value.
#' @name integer
#'
#' @examples
#' vec <- as.integer(fm.runif(1000, min=0, max=10))
NULL

#' @rdname integer
setMethod("as.integer", "fm", function(x) {
		  if (.typeof.int(x) == "integer")
			  x
		  else
			  fm.sapply(x, fm.buo.as.int)
	})
#' @rdname integer
setMethod("as.integer", "fmV", function(x) {
		  if (.typeof.int(x) == "integer")
			  x
		  else
			  fm.sapply(x, fm.buo.as.int)
	})

#' Numeric Vectors
#'
#' Coerces or test objects of type \code{"numeric"}. \code{is.numeric}
#' is a more general test of an object being interpretable as numbers.
#'
#' @param x a FlashR object to be coerced or tested.
#' @return \code{as.numeric} returns a numeric FlashR object,
#' \code{is.numeric} returns a logical value.
#' @name numeric
#'
#' @examples
#' vec <- as.numeric(fm.runif(1000, min=0, max=10))
#' is.numeric(vec)
NULL

#' @rdname numeric
setMethod("as.numeric", "fm", function(x) {
		  if (.typeof.int(x) == "double")
			  x
		  else
			  fm.sapply(x, fm.buo.as.numeric)
	})
#' @rdname numeric
setMethod("as.numeric", "fmV", function(x) {
		  if (.typeof.int(x) == "double")
			  x
		  else
			  fm.sapply(x, fm.buo.as.numeric)
	})
#' @rdname numeric
setMethod("is.numeric", "fm", function(x)
		  .typeof.int(x) == "double" || .typeof.int(x) == "integer")
#' @rdname numeric
setMethod("is.numeric", "fmV", function(x)
		  .typeof.int(x) == "double" || .typeof.int(x) == "integer")

.fmV2scalar <- function(x)
{
	fm.conv.FM2R(x)[1]
}

#' Extract Parts of an object
#'
#' Operators acting on FlashR vectors and matrices to extract parts of
#' the objects.
#'
#' FlashR supports extracting data from an object as well as replacing data
#' of a matrix. To assign data in a matrix, it only supports setting
#' the entire rows or columns.
#'
#' @param x object from which to extract elements.
#' @param i,j indices specifying elements to extract. Indices are integer
#'        or numeric vectors.
#' @param drop for matrices. If \code{TRUE}, the result is coerced to a vector.
#' @param value a FlashR vector or matrix that contains new data for the rows
#"        or columns.
#' @name Extract
#'
#' @examples
#' mat <- fm.runif.matrix(100, 10)
#' mat[1:5,]
#' mat[,1:5]
#' mat[1:5,1:5]
#' mat[,1] <- fm.runif(nrow(mat))
NULL

#' @rdname Extract
setMethod("[", signature(x="fm", j="missing"),
		  function(x, i, j, drop=TRUE) {
			  if (is.matrix(i) || fm.is.matrix(i)) {
				  print("doesn't support boolean index matrix")
				  NULL
			  }
			  else {
				  ret <- fm.get.rows(x, i)
				  if (drop && length(i) == 1)
					  ret <- fm.as.vector(ret)
				  ret
			  }
		  })
#' @rdname Extract
setMethod("[", signature(x="fm", i="missing"),
		  function(x, i, j, drop=TRUE) {
			  ret <- fm.get.cols(x, j)
			  if (drop && length(j) == 1)
				  ret <- fm.as.vector(ret)
			  ret
		  })
#' @rdname Extract
setMethod("[", signature(x="fm"), function(x, i, j, drop=TRUE) {
		  if (nrow(x) > ncol(x)) {
			  ret <- fm.get.cols(x, j)
			  ret <- fm.get.rows(ret, i)
		  }
		  else {
			  ret <- fm.get.rows(x, i)
			  ret <- fm.get.cols(ret, j)
		  }
		  if (drop && (length(i) == 1 || length(j) == 1))
			  ret <- fm.as.vector(ret)
		  ret
		  })
#' @rdname Extract
setMethod("[", signature(x="fmV"), function(x, i) fm.get.eles.vec(x, i))
#' @rdname Extract
setMethod("[<-", signature(x="fm", j="missing"),
		  function(x, i, j, value) fm.set.rows(x, i, value))
#' @rdname Extract
setMethod("[<-", signature(x="fm", i="missing"),
		  function(x, i, j, value) fm.set.cols(x, j, value))

#' Return the First or Last part of an Object
#'
#' Returns the first or last parts of a FlashR vector or matrix.
#'
#' @param x a FlashR vector or matrix.
#' @param n a positive integer.
#' @return An object (usually) like \code{x} but generally smaller.
#' @name head
#'
#' @examples
#' mat <- fm.runif.matrix(100, 10)
#' head(mat)
#' tail(mat)
NULL

#' @rdname head
setMethod("head", signature(x="fm"), function(x, n=6L) fm.get.rows(x, 1:n))
#' @rdname head
setMethod("tail", signature(x="fm"), function(x, n=6L) fm.get.rows(x, (nrow(x)-n+1):nrow(x)))
#' @rdname head
setMethod("head", signature(x="fmV"), function(x, n=6L) x[1:n])
#' @rdname head
setMethod("tail", signature(x="fmV"), function(x, n=6L) x[(length(x)-n+1):length(x)])

#' Sweep out Array Summaries
#'
#' Return an array obtained from an input array by sweeping out a
#' summary statistic.
#' @param x a FlashR matrix.
#' @param MARGIN an integer giving the extent of \code{x} which
#'        correspond to \code{STATS}.
#' @param STATS the summary statistic which is to be swept out.
#' @param FUN the function to be used to carry out the sweep.
#' @param check.margin logical. It's ignored right now.
#' @param ... optional arguments to \code{FUN}.
#' @return A matrix with the same shape as \code{x}, but with the summary
#' statistics swept out.
#' @name sweep
#'
#' @examples
#' mat <- fm.runif.matrix(100, 10)
#' sweep(mat, 1, runif(100))
setMethod("sweep", "fm",
		  function(x, MARGIN, STATS, FUN="-", check.margin=TRUE, ...) {
			  if (MARGIN == 2)
				  fm.mapply.row(x, STATS, FUN)
			  else if (MARGIN == 1)
				  fm.mapply.col(x, STATS, FUN)
			  else
				  stop("a wrong margin")
		  })

#' Matrix Crossproduct
#'
#' Given matrices \code{x} and \code{y} as arguments, return a matrix
#' cross-product. This is formally equivalent to (but usually slightly
#' faster than) the call \code{t(x) \%*\% y} (\code{crossprod}) or
#' \code{x \%*\% t(y)} (\code{tcrossprod}).
#'
#' @param x,y numeric matrices (or vectors): \code{y = NULL} is taken
#' to be the same matrix as \code{x}. Vectors are promoted to single-column
#' or single-row matrices, depending on the context.
#' @param lazy a logical value, indicating to perform the computation lazily.
#' @return a double matrix
#' @name crossprod
#'
#' @examples
#' mat <- fm.runif.matrix(100, 10)
#' crossprod(mat)
#' tcrossprod(mat)
NULL

#' @rdname crossprod
setMethod("crossprod", "fm", function(x, y=NULL) {
		  if (is.null(y))
			  y <- x
		  fm.multiply(t(x), y)
})
#' @rdname crossprod
setMethod("tcrossprod", "fm", function(x, y=NULL) {
		  if (is.null(y))
			  y <- x
		  fm.multiply(x, t(y))
})

#' FlashR Summaries
#'
#' \code{summary} produces summaries of a FlashR object.
#'
#' It computes the min, max, sum, mean, L1, L2, number of non-zero values if
#' the argument is a vector, or these statistics for each column if the argument
#' is a matrix.
#'
#' @param x a FlashR vector or matrix.
#' @return A list containing the following named components:
#' \itemize{
#' \item{min}{The minimum value}
#' \item{max}{The maximum value}
#' \item{mean}{The average}
#' \item{normL1}{The L1 norm}
#' \item{normL2}{The L2 norm}
#' \item{numNonzeros}{The number of non-zero values}
#' }
#' @name summary
#'
#' @examples
#' mat <- fm.runif.matrix(100, 10)
#' summary(mat)
NULL

#' @rdname summary
.summary <- function(x)
{
	x <- fm.as.matrix(x)
	res <- list()
	res[[1]] <- fm.agg.mat(x, 2, fm.bo.min)
	res[[2]] <- fm.agg.mat(x, 2, fm.bo.max)
	res[[3]] <- fm.agg.mat(x, 2, fm.bo.add)
	res[[4]] <- fm.agg.mat(abs(x), 2, fm.bo.add)
	res[[5]] <- fm.agg.mat(x * x, 2, fm.bo.add)
	res[[6]] <- fm.agg.mat(x != 0, 2, fm.bo.add)
	res <- lapply(res, function(o) fm.conv.FM2R(o))
	mean <- res[[3]]/nrow(x)
	var <- (res[[5]]/nrow(x) - mean^2) * nrow(x) / (nrow(x) - 1)
	list(min=res[[1]], max=res[[2]], mean=mean, normL1=res[[4]],
		 normL2=sqrt(res[[5]]), numNonzeros=res[[6]], var=var)
}

#' @rdname summary
setMethod("summary", "fm", function(object, ...) .summary(object))
#' @rdname summary
setMethod("summary", "fmV", function(object, ...) .summary(object))

#' Eigensolver
#'
#' \code{fm.eigen} computes eigenvalues/vectors of a square matrix.
#' \code{fm.cal.residul} computes the residual of the eigenvalues.
#'
#' \code{fm.eigen} uses Anasazi package of Trilinos, if Anasazi is compiled
#' into FlashR, or eigs to compute eigenvalues.
#'
#' The \code{which} specify which eigenvalues/vectors to compute, character
#' constant with exactly two characters. Possible values for symmetric input
#' matrices:
#' \itemize{
#' \item{"LR"}{Compute \code{k} largest (real part) eigenvalues.}
#' \item{"SR"}{Compute \code{k} smallest (real part) eigenvalues.}
#' \item{"LM"}{Compute \code{k} largest (in magnitude) eigenvalues.}
#' \item{"SM"}{Compute \code{k} smallest (in magnitude) eigenvalues.}
#' }
#'
#' @param mul The function to perform the matrix-vector multiplication.
#' @param k Integer. The number of eigenvalues to compute.
#' @param which String. Selection criteria.
#' @param sym Logical scalar, whether the input matrix is symmetric.
#' @param options List. Additional options to the eigensolver.
#' @param env The environment in which \code{mul} will bevaluated.
#' @param values The eigenvalues
#' @param vectors The eigenvectors
#'
#' The \code{options} argument specifies what kind of computation to perform.
#' It is a list with the following members, which correspond directly to
#' Anasazi parameters:
#'
#' solver String. The name of the eigensolver to solve the eigenproblems.
#'				Currently, it supports three eigensolvers: KrylovSchur,
#'				Davidson and LOBPCG. KrylovSchur is the default eigensolver.
#'
#' tol Numeric scalar. Stopping criterion: the relative accuracy of
#'				the Ritz value is considered acceptable if its error is less
#'				than \code{tol} times its estimated value.
#'
#' block_size Numeric scalar. The eigensolvers use a block extension of an
#'				eigensolver algorithm. The block size determines the number
#'				of the vectors that operate together.
#'
#' num_blocks Numeric scalar. The number of blocks to compute eigenpairs.
#'
#'
#' @return \code{fm.eigen} returns a named list with the following members:
#'         values: Numeric vector, the desired eigenvalues.
#'         vectors: Numeric matrix, the desired eigenvectors as columns.
#'         \code{fm.cal.residul} returns the corresponding residuals for
#'         the eigenvalues.
#' @name fm.eigen
#' @author Da Zheng <dzheng5@@jhu.edu>
#'
#' @examples
#' mat <- fm.load.sparse.matrix("./spm123.mat", "./spm123.mat_idx")
#" mul <- function(x, extra) mat %*% x
#' res <- fm.eigen(mul, 10, nrow(mat))
#' fm.cal.residul(mul, res$values, res$vectors)
fm.eigen <- function(mul, k, n, which=c("LM", "SM", "LR", "SR"),
					 sym=TRUE, options=NULL, env = parent.frame())
{
	if (!sym) {
		print("fm.eigen only supports symmetric matrices")
		return(NULL)
	}
	if (k > n) {
		print("k should be smaller than n")
		return(NULL)
	}
	which <- match.arg(which)
	# We only want to solve a very large eigenvalue problem with Anasazi.
	if (is.loaded("R_FM_eigen") && n > 1000000) {
		if (is.null(options))
			options <- list()
		options$nev <- k
		options$n <- n
		options$which <- which
		ret <- .Call("R_FM_eigen", as.function(mul), NULL, as.logical(sym),
					 as.list(options), PACKAGE="FlashR")
		ret$vectors <- .new.fm(ret$vectors)
		ret
	}
	else {
		eigs.opts <- list()
		if (!is.null(options)) {
			if (!is.null(options[["tol"]]))
				eigs.opts$tol <- options$tol
			if (!is.null(options[["block_size"]])
				&& !is.null(options[["num_blocks"]]))
				eigs.opts$ncv <- options$block_size * options$num_blocks
			# If we compute a small number of eigenvalues, we need a larger
			# subspace.
			else if (k <= 2)
				eigs.opts$ncv <- k * 2 + 3
			else
				eigs.opts$ncv <- k * 2
		}
		eigs.fun <- function(x, extra) {
			# TODO this might be slow
			fm.x <- fm.as.vector(x)
			ret <- mul(fm.x, extra)
			# TODO this might be slow
			as.vector(ret)
		}
		if (sym)
			eigs.ret <- eigs_sym(eigs.fun, k=k, which=which, n=n, opts=eigs.opts)
		else
			eigs.ret <- eigs(eigs.fun, k=k, which=which, n=n, opts=eigs.opts)
		list(values=eigs.ret$values, vectors=fm.as.matrix(eigs.ret$vectors))
	}
}

#' @rdname fm.eigen
fm.cal.residul <- function(mul, values, vectors)
{
	tmp <- mul(vectors, NULL) - fm.mapply.row(vectors, values, "*")
	l2 <- sqrt(colSums(tmp * tmp))
	fm.conv.FM2R(l2) / values
}

#' Scaling and Centering of Matrix
#'
#' \code{scale} centers and/or scales the columns of a FlashR matrix.
#'
#' The value of \code{center} determines how column centering is
#' performed.  If \code{center} is a numeric vector with length equal to
#' the number of columns of \code{x}, then each column of \code{x} has the
#' corresponding value from \code{center} subtracted from it.  If \code{center}
#' is \code{TRUE} then centering is done by subtracting the column means
#' (omitting \code{NA}s) of \code{x} from their corresponding columns, and if
#' \code{center} is \code{FALSE}, no centering is done.
#'
#' The value of \code{scale} determines how column scaling is performed
#' (after centering).  If \code{scale} is a numeric vector with length
#' equal to the number of columns of \code{x}, then each column of \code{x} is
#' divided by the corresponding value from \code{scale}.  If \code{scale} is
#' \code{TRUE} then scaling is done by dividing the (centered) columns of
#' \code{x} by their standard deviations if \code{center} is \code{TRUE},
#' and the root mean square otherwise.  If \code{scale} is \code{FALSE},
#' no scaling is done.
#'
#' The root-mean-square for a (possibly centered) column is defined
#' as sqrt(sum(x^2)/(n-1)), where x is a vector of the non-missing
#' values and n is the number of non-missing values.  In the case
#' \code{center = TRUE}, this is the same as the standard deviation, but
#' in general it is not.
#'
#' @param x a FlashR matrix
#' @param center either a logical value or a numeric vector of length equal to
#'        the number of columns of \code{x}.
#' @param sclae either a logical value or a numeric vector of length equal to
#'        the number of columns of \code{x}.
#' @return a FlashR matrix.
#' @author Da Zheng <dzheng5@@jhu.edu>
#' @name scale
#'
#' @examples
#' mat <- fm.runif.matrix(100, 10)
#' mat <- scale(mat)
setMethod("scale", "fm", function(x, center=TRUE, scale=TRUE) {
		  # The output of this function is floating-point values.
		  # It's better to convert the element type in advance to avoid
		  # any overflow.
		  x <- as.double(x)
		  # TODO it needs to handle NA.
		  # If the center is true, center columns by their means.
		  if (!is.logical(center)) {
			  if (length(center) != ncol(x)) {
				  print("The length of center should be equal to #columns of x")
				  return(NULL)
			  }
			  x <- fm.mapply.row(x, center, fm.bo.sub)
		  }
		  if (!is.logical(scale)) {
			  if (length(scale) != ncol(x)) {
				  print("The length of scale should be equal to #columns of x")
				  return(NULL)
			  }
			  x <- fm.mapply.row(x, scale, fm.bo.div)
		  }

		  # If scale is true and center is also true, scale by their standard
		  # deviation.
		  if (is.logical(scale) && scale && is.logical(center) && center) {
			  sum.x <- colSums(x)
			  sum.x2 <- colSums(x * x)
			  n <- nrow(x)
			  avg <- sum.x / n
			  sd <- sqrt((sum.x2 - n * avg * avg) / (n - 1))
			  center <- sum.x / n
			  x <- fm.mapply.row(x, center, fm.bo.sub)
			  x <- fm.mapply.row(x, sd, fm.bo.div)
		  }
		  # If scale is true and center is false, scale by their root mean
		  # square.
		  else if (is.logical(scale) && scale) {
			  sum.x2 <- colSums(x * x)
			  scal <- sqrt(sum.x2 / (nrow(x) - 1))
			  x <- fm.mapply.row(x, scal, fm.bo.div)
		  }
		  # If scale is false and center is true.
		  else if (is.logical(center) && center) {
			  center <- colMeans(x)
			  x <- fm.mapply.row(x, center, fm.bo.sub)
		  }

		  x
})

#' Drop Redundant Extent Information
#'
#' Delete the dimensions of a FlashR matrix which have only one level.
#'
#' If the input matrix has only one row or one column, it works the same as
#' fm.as.vector. Otherwise, it returns the original matrix.
#'
#' @param x a FlashR matrix.
#' @return a FlashR matrix or vector.
#' @name drop
#'
#' @examples
#' mat <- fm.runif.matrix(100, 1)
#' vec <- drop(mat)
setMethod("drop", "fm", function(x) {
		  if (nrow(x) > 1 && ncol(x) > 1)
			  x
		  else
			  fm.as.vector(x)
})

#' Which indices are TRUE?
#'
#' Give the \code{TRUE} indices of a logical object, allowing for array
#' indices.
#'
#' @param x a \code{logical} vector or array. \code{NA}s are allowed and omitted
#'        (treated as if \code{FALSE}).
#' @param arr.ind logical; should *arr*ay *ind*ices be returned when \code{x} is
#'        an array. This argument is ignored currently.
#' @param useNames logical indicating if the value of \code{arrayInd()} should
#'        have (non-null) dimnames at all. This argument is ignored currently.
#' @return a FlashR vector
#' @name which
#'
#' @examples
#' vec <- fm.runif(1000000)
#' idx <- which(vec < 0.5)
NULL

#' @rdname which
setMethod("which", "fm", function(x, arr.ind=FALSE, useNames=TRUE) {
		  fm.seq.int(1, length(x), 1)[as.logical(x)]
})
#' @rdname which
setMethod("which", "fmV", function(x, arr.ind=FALSE, useNames=TRUE) {
		  fm.seq.int(1, length(x), 1)[as.logical(x)]
})
