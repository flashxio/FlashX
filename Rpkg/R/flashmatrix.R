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

setClass("fm", representation(pointer = "externalptr", name = "character",
							  nrow = "numeric", ncol = "numeric",
							  type="character", ele_type="character"))
setClass("fmV", representation(pointer = "externalptr", name = "character",
							   len = "numeric", type="character",
							   ele_type="character"))
setClass("fmFactorV", representation(num.levels = "integer"), contains = "fmV")
# This represents a lazily evaluated agg matrix, wide inner product matrix and
# groupby matrix.
setClass("fmSink", representation(pointer = "externalptr", name = "character",
								  nrow = "numeric", ncol = "numeric",
								  type="character", ele_type="character"))
# We use symbolic representation for UDFs.
# We will select the right form and right element type for a UDF when
# the UDF is used.
setClass("fm.bo", representation(info = "integer", name = "character"))
setClass("fm.agg.op", representation(agg = "integer", combine = "integer",
									 name = "character"))

new.fm <- function(fm)
{
	if (!is.null(fm))
		new("fm", pointer=fm$pointer, name=fm$name, nrow=fm$nrow, ncol=fm$ncol,
			type=fm$type, ele_type=fm$ele_type)
	else
		NULL
}

new.fmV <- function(fm)
{
	if (!is.null(fm))
		new("fmV", pointer=fm$pointer, name=fm$name, len=fm$len, type=fm$type,
			ele_type=fm$ele_type)
	else
		NULL
}

new.fmFactorV <- function(fm)
{
	if (!is.null(fm))
		new("fmFactorV", num.levels=fm$levels, pointer=fm$pointer,
			name=fm$name, len=fm$len, type=fm$type, ele_type=fm$ele_type)
	else
		NULL
}

new.fmSink <- function(fm)
{
	if (!is.null(fm))
		new("fmSink", pointer=fm$pointer, name=fm$name, nrow=fm$nrow, ncol=fm$ncol,
			type=fm$type, ele_type=fm$ele_type)
}

#' Reconfigure FlashMatrix
#'
#' This reconfigures FlashMatrix with the settings in the configuration file.
#' The configuration file contains a list of key-value pairs. Each line in
#' the file is a key-value pair in the form of "key_name=value".
#' @param conf.file The configuration file.
#' @name fm.set.conf
#' @author Da Zheng <dzheng5@@jhu.edu>
fm.set.conf <- function(conf.file)
{
	ret <- .Call("R_FM_set_conf", conf.file, PACKAGE="FlashR")
	stopifnot(ret);
}

fm.set.log.level <- function(level)
{
	.Call("R_FM_set_log_level", level, PACKAGE="FlashR")
}

fm.print.conf <- function()
{
	fg.set.log.level("info")
	.Call("R_SAFS_print_conf", PACKAGE="FlashR")
	.Call("R_FM_print_conf", PACKAGE="FlashR")
	fg.set.log.level("warning")
}

#' Indicate whether a matrix has been loaded to FlashR
#'
#' This function indicates whether a matrix has been loaded to FlashR.
#' @param name A matrix name.
#' @return A boolean value: true if the matrix has been loaded to FlashR.
#' @name fm.exist.matrix
#' @author Da Zheng <dzheng5@@jhu.edu>
fm.exist.matrix <- function(name)
{
	.Call("R_FM_exist_matrix", name, PACKAGE="FlashR")
}

fm.get.params <- function(name)
{
	.Call("R_FM_get_params", name, PACKAGE="FlashR")
}

#' Load a sparse matrix to FlashR.
#'
#' Load a sparse matrix to FlashR from difference sources.
#' 
#' `fm.get.sparse.matrix' gets a FlashMatrix sparse matrix that references
#' a matrix represented by a FlashGraph object.
#'
#' `fm.get.dense.matrix' gets a FlashMatrix dense matrix from the underlying
#" filesystem.
#'
#' `fm.load.sparse.matrix' loads a FlashMatrix sparse matrix from files.
#' The matrix in the file is in the FlashMatrix format.
#'
#' @param fg A FlashGraph object.
#' @param mat.file The file that stores the sparse matrix.
#' @param index.file The file that stores the index of the sparse matrix.
#' @return a FlashMatrix matrix.
#' @name fm.get.matrix
#' @author Da Zheng <dzheng5@@jhu.edu>
#' @examples
#' fg <- fg.load.graph("graph.adj", "graph.index")
#' fm <- fm.get.sparse.matrix(fg)
fm.get.sparse.matrix <- function(fg)
{
	stopifnot(!is.null(fg))
	stopifnot(class(fg) == "fg")
	stopifnot(fg.exist.graph(fg$name))
	m <- .Call("R_FM_get_matrix_fg", fg, PACKAGE="FlashR")
	new.fm(m)
}

#' @rdname fm.get.matrix
fm.get.dense.matrix <- function(name)
{
	stopifnot(!is.null(name))
	stopifnot(class(name) == "character")
	m <- .Call("R_FM_get_dense_matrix", name, PACKAGE="FlashR")
	new.fm(m)
}

#' @rdname fm.get.matrix
fm.load.sparse.matrix <- function(mat, index, t.mat=NULL, t.index=NULL, in.mem=TRUE)
{
	if (is.null(t.mat) || is.null(t.index))
		m <- .Call("R_FM_load_matrix_sym", mat, index, in.mem, PACKAGE="FlashR")
	else
		m <- .Call("R_FM_load_matrix_asym", mat, index, t.mat, t.index, in.mem,
			  PACKAGE="FlashR")
	new.fm(m)
}

#' Create a FlashMatrix vector with replicated elements.
#'
#' @param x A constant initial value
#' @param times The length of the vector to be generated.
#' @return A FlashMatrix vector
#' @name fm.rep.int
#' @author Da Zheng <dzheng5@@jhu.edu>
fm.rep.int <- function(x, times)
{
	if (times <= 0) {
		print("we can't generate a vector of 0 elements")
		return(NULL)
	}
	stopifnot(is.vector(x) && is.atomic(x))
	vec <- .Call("R_FM_create_vector", as.numeric(times), x, PACKAGE="FlashR")
	new.fmV(vec)
}

#' Create a FlashMatrix vector with a sequence of numbers.
#'
#' @param from the starting value of the sequence.
#' @param to the end value of the sequence.
#' @param by number: increment of the sequence.
#' @return a sequence of numbers that belongs to [from, to]
#' @name fm.seq.int
#' @author Da Zheng <dzheng5@@jhu.edu>
fm.seq.int <- function(from, to, by)
{
	if ((to - from) / by <= 0) {
		print("we can't generate a vector of 0 elements")
		return(NULL)
	}
	vec <- .Call("R_FM_create_seq", from, to, by, PACKAGE="FlashR")
	new.fmV(vec)
}

#' Create a FlashMatrix vector with uniformly random numbers.
#'
#' @param n the number of random numbers to be generated.
#' @param min lower limits of the distribution.
#' @param max upper limits of the distribution.
#' @name fm.runif
#' @author Da Zheng <dzheng5@@jhu.edu>
fm.runif <- function(n, min=0, max=1)
{
	if (n <= 0) {
		print("we can't generate a vector of 0 elements")
		return(NULL)
	}
	vec <- .Call("R_FM_create_rand", n, min, max, PACKAGE="FlashR")
	new.fmV(vec)
}

setMethod("as.vector", signature(x = "fmV"), function(x) fm.conv.FM2R(x))
setMethod("as.matrix", signature(x = "fm"), function(x) fm.conv.FM2R(x))

#' Convert the data layout of a FlashMatrix matrix.
#'
#' A FMR matrix can store elements in a row-major or column-major order.
#' By changing the data layout, we can improve efficiency of some matrix
#' operations.
#'
#' @param fm a FlashMatrix matrix
#' @param byrow a logical value to determine the data layout of a FlashMatrix
#' matrix.
#' @return a FlashMatrix matrix
#' @author Da Zheng <dzheng5@@jhu.edu>
fm.conv.layout <- function(fm, byrow=FALSE)
{
	stopifnot(!is.null(fm))
	stopifnot(class(fm) == "fm")
	ret <- .Call("R_FM_conv_layout", fm, as.logical(byrow), PACKAGE="FlashR")
	if (!is.null(ret))
		new.fm(ret)
	else
		NULL
}

#' Get the data layout of a FlashMatrix matrix.
#'
#' @param fm a FlashMatrix matrix
#' @return a string that indicates the data layout.
#' @author Da Zheng <dzheng5@@jhu.edu>
fm.get.layout <- function(fm)
{
	stopifnot(!is.null(fm))
	stopifnot(class(fm) == "fm")
	.Call("R_FM_get_layout", fm, PACKAGE="FlashR")
}

#' Convert a regular R object to a FlashMatrix object.
#'
#' If the R object is a matrix, `byrow' determines how data in the generated
#' FlashMatrix object is organized in memory.
#'
#' @param obj a regular R object
#' @param byrow a logical value to determine the data layout of a FlashMatrix
#' matrix.
#' @return a FlashMatrix object.
#' @name fm.conv.R2FM
#' @author Da Zheng <dzheng5@@jhu.edu>
fm.conv.R2FM <- function(obj, byrow=FALSE)
{
	stopifnot(!is.null(obj))
	# This function only deals with vectors and matrices of primitive types.
	stopifnot(is.atomic(obj))

	if (is.vector(obj)) {
		vec <- .Call("R_FM_conv_RVec2FM", obj, PACKAGE="FlashR")
		new.fmV(vec)
	}
	else {
		m <- .Call("R_FM_conv_RMat2FM", obj, as.logical(byrow), PACKAGE="FlashR")
		new.fm(m)
	}
}

#' Convert a FlashMatrix object to a regular R object
#'
#' @param obj a FlashMatrix object
#' @return a regular R object.
#' @name fm.conv.FM2R
#' @author Da Zheng <dzheng5@@jhu.edu>
fm.conv.FM2R <- function(obj)
{
	stopifnot(!is.null(obj))
	if (class(obj) == "fm") {
		nrow <- dim(obj)[1]
		ncol <- dim(obj)[2]
		if (fm.typeof(obj) == "integer")
			ret <- matrix(vector(mode="integer", nrow * ncol), nrow, ncol)
		else if (fm.typeof(obj) == "double")
			ret <- matrix(vector(mode="double", nrow * ncol), nrow, ncol)
		else if (fm.typeof(obj) == "logical")
			ret <- matrix(vector(mode="logical", nrow * ncol), nrow, ncol)
		.Call("R_FM_copy_FM2R", obj, ret, PACKAGE="FlashR")
		ret
	}
	else if (class(obj) == "fmV") {
		len <- length(obj)
		if (fm.typeof(obj) == "integer")
			ret <- vector(mode="integer", len)
		else if (fm.typeof(obj) == "double")
			ret <- vector(mode="double", len)
		else if (fm.typeof(obj) == "logical")
			ret <- vector(mode="logical", len)
		.Call("R_FM_copy_FM2R", obj, ret, PACKAGE="FlashR")
		ret
	}
	else {
		print("It has to be a FlashMatrix object")
		NULL
	}
}

#' Create a matrix from the given set of values.
#'
#' The given set has to be a FlashMatrix vector.
#'
#' @param vec a given vector.
#' @param nrow the number of rows in the output matrix.
#' @param ncol the number of columns in the output matrix.
#' @param byrow logical. If FALSE (the default) the matrix is filly
#'				columns, otherwise the matrix is filled by rows.
#' @return a FlashMatrix matrix
fm.matrix <- function(vec, nrow, ncol, byrow=FALSE)
{
	stopifnot(!is.null(vec))
	stopifnot(class(vec) == "fmV")
	m <- .Call("R_FM_conv_matrix", vec, as.numeric(nrow), as.numeric(ncol),
				 as.logical(byrow), PACKAGE="FlashR")
	new.fm(m)
}

#' The information of a FlashMatrix object
#'
#' Functions for providing the basic information of a matrix.
#'
#' `fm.is.sym' indicates whether a matrix is symmetric.
#'
#' `fm.matrix.layout' indicates how data in a matrix is organized.
#
#' `fm.is.sparse' indicates whether a matrix is sparse.
#'
#' `fm.is.vector' indicates whether a FlashMatrix object is a vector.
#' `fm.is.matrix' indicates whether a FlashMatrix object is a matrix.
#'
#' @param fm The FlashMatrix object
#' @return `fm.is.sym' and `fm.is.sparse' returns boolean constants.
#' @name fm.info
#' @author Da Zheng <dzheng5@@jhu.edu>

#' @rdname fm.info
fm.is.sym <- function(fm)
{
	stopifnot(!is.null(fm))
	stopifnot(class(fm) == "fm")
	.Call("R_FM_is_sym", fm, PACKAGE="FlashR")
}

#' @rdname fm.info
fm.matrix.layout <- function(fm)
{
	stopifnot(!is.null(fm))
	stopifnot(class(fm) == "fm")
	.Call("R_FM_matrix_layout", fm, PACKAGE="FlashR")
}

#' @rdname fm.info
fm.is.sparse <- function(fm)
{
	stopifnot(!is.null(fm))
	stopifnot(class(fm) == "fm")
	fm@type == "sparse"
}

#' @rdname fm.info
fm.is.vector <- function(fm)
{
	stopifnot(!is.null(fm))
	class(fm) == "fmV"
}

#' @rdname fm.info
fm.is.matrix <- function(fm)
{
	stopifnot(!is.null(fm))
	class(fm) == "fm"
}

fm.typeof <- function(fm)
{
	stopifnot(!is.null(fm))
	stopifnot(class(fm) == "fm" || class(fm) == "fmV")
	fm@ele_type
}

#' Convert a FlashMatrix matrix to a FlashMatrix vector.
#'
#' The matrix must have only one row or one column. Otherwise, the function
#' returns an error.
#'
#' @param fm a FlashMatrix matrix
#' @return a FlashMatrix vector
fm.as.vector <- function(obj)
{
	stopifnot(!is.null(obj))
	if (class(obj) == "fmV")
		vec
	else if (class(obj) == "fm") {
		vec <- .Call("R_FM_as_vector", obj, PACKAGE="FlashR")
		if (!is.null(vec))
			new.fmV(vec)
		else
			NULL
	}
	else if (is.vector(obj))
		fm.conv.R2FM(obj)
	else
		NULL
}

fm.as.matrix <- function(obj)
{
	stopifnot(!is.null(obj))
	if (class(obj) == "fm")
		obj
	else if (class(obj) == "fmV") {
		# A FlashMatrix vector is actually stored in a dense matrix.
		# We only need to construct the fm object in R.
		new("fm", pointer=obj@pointer, name=obj@name, nrow=obj@len,
			ncol=1, type=obj@type)
	}
	else {
		# Let's convert it to a FM object
		ret <- fm.conv.R2FM(obj)
		# Then try to convert it to a FM matrix.
		fm.as.matrix(ret)
	}
}

#' Convert a FlashMatrix vector to a FlashMatrix factor vector.
#'
#' @param fm a FlashMatrix vector.
#' @return a FlashMatrix factor vector.
fm.as.factor <- function(fm, num.levels = -1)
{
	stopifnot(!is.null(fm))
	stopifnot(class(fm) == "fmV")
	vec <- .Call("R_FM_as_factor_vector", fm, as.integer(num.levels),
				 PACKAGE="FlashR")
	new.fmFactorV(vec)
}

#' Matrix multiplication
#'
#' Multiply a sparse/dense matrix with a dense vector/matrix.
#'
#' @param fm A FlashMatrix matrix
#' @param mat A FlashMatrix dense matrix.
#' @return a FlashMatrix vector if the second argument is a vector;
#' a FlashMatrix matrix if the second argument is a matrix.
#' @name fm.multiply
#' @author Da Zheng <dzheng5@@jhu.edu>
fm.multiply <- function(fm, mat)
{
	stopifnot(!is.null(fm) && !is.null(mat))
	stopifnot(class(fm) == "fm")
	if (class(mat) == "fmV") {
		stopifnot(dim(fm)[2] == length(mat))
	}
	else {
		stopifnot(!fm.is.sparse(mat))
		stopifnot(dim(fm)[2] == dim(mat)[1])
	}

	if (fm.is.sparse(fm))
		o <- .Call("R_FM_multiply_sparse", fm, mat, PACKAGE="FlashR")
	else
		o <- .Call("R_FM_multiply_dense", fm, mat, PACKAGE="FlashR")
	if (class(mat) == "fmV")
		new.fmV(o)
	else
		new.fm(o)
}

#' Matrix inner product
#'
#' It takes two operators and performs inner product on a dense matrix
#" and a dense vector/matrix.
#'
#' @param fm A FlashMatrix matrix
#' @param mat A FlashMatrix dense matrix.
#' @param Fun1 The reference or the name of one of the predefined basic binary
#' operators.
#' @param Fun2 The reference or the name of one of the predefined basic binary
#' operators.
#' @return a FlashMatrix vector if the second argument is a vector;
#' a FlashMatrix matrix if the second argument is a matrix.
#' @name fm.inner.prod
#' @author Da Zheng <dzheng5@@jhu.edu>
fm.inner.prod <- function(fm, mat, Fun1, Fun2)
{
	stopifnot(!is.null(fm) && !is.null(mat))
	stopifnot(class(fm) == "fm")
	if (class(mat) == "fmV") {
		stopifnot(dim(fm)[2] == length(mat))
	}
	else {
		stopifnot(!fm.is.sparse(mat))
		stopifnot(dim(fm)[2] == dim(mat)[1])
	}
	if (class(Fun1) == "character")
		Fun1 <- fm.get.basic.op(Fun1)
	stopifnot(class(Fun1) == "fm.bo")
	if (class(Fun2) == "character")
		Fun2 <- fm.get.basic.op(Fun2)
	stopifnot(class(Fun2) == "fm.bo")

	if (fm.is.sparse(fm)) {
		print("inner product doesn't support sparse matrices yet")
		return(NULL)
	}
	else
		o <- .Call("R_FM_inner_prod_dense", fm, mat, Fun1, Fun2,
				   PACKAGE="FlashR")
	if (class(mat) == "fmV")
		new.fmV(o)
	else
		new.fm(o)
}

#' The basic operators supported by FlashMatrix.
#'
#'
#' The basic operators are mainly used by the FlashMatrix functions that
#' accept operators as arguments. Such a function includes `fm.mapply',
#' `fm.inner.prod', etc.
#'
#' `fm.get.basic.op' gets the predefined basic binary operator specified by a user.
#' The supported basic binary operators are:
#' \itemize{
#' \item{"add" or "+"}{compute addition.}
#' \item{"sub" or "-"}{compute subtraction;}
#' \item{"mul" or "*"}{compute multiplication;}
#' \item{"div" or "/"}{compute division;}
#' \item{"min" and "max"}{compute minimum and maximum, respectively;}
#' \item{"pow"}{compute exponential;}
#' \item{"eq" or "=="}{compute equality;}
#' \item{"gt" or ">"}{compute greater than;}
#' \item{"ge" or ">="}{compute greater than or equal to;}
#' \item{"lt" or "<"}{compute less than;}
#' \item{"le" or "<="}{compute less than or equal to;}
#' }
#'
#' `fm.get.basic.uop' gets the predefined basic unary operator specified by a user.
#' The supported basic unary operators are:
#' \itemize{
#' \item{"neg"}{compute negate;}
#' \item{"sqrt"}{compute square root;}
#' \item{"abs"}{compute absolute value;}
#' \item{"not"}{compute logical NOT;}
#' \item{"ceil" and "floor"}{compute a ceiling and a floor, respectively;}
#' \item{"log", "log2" and "log10"}{compute log with different bases;}
#' \item{"round"}{round a number;}
#' \item{"as.int" and "as.numeric"}{cast a number to an integer and a numeric
#' value, respectively.}
#' }
#'
#' `fm.init.basic.op' initializes the following basic operators.
#' \itemize{
#' \item{`fm.bo.add'}{the predifined basic binary operator for addition.}
#' \item{`fm.bo.sub'}{the predifined basic binary operator for subtraction.}
#' \item{`fm.bo.mul'}{the predifined basic binary operator for multiplication.}
#' \item{`fm.bo.div'}{the predifined basic binary operator for division.}
#' \item{`fm.bo.min'}{the predifined basic binary operator for computing minimum.}
#' \item{`fm.bo.max'}{the predifined basic binary operator for computing maximum.}
#' \item{`fm.bo.pow'}{the predifined basic binary operator for computing exponential.}
#' \item{`fm.bo.eq', `fm.bo.neq', `fm.bo.gt', `fm.bo.ge', `fm.bo.lt' and `fm.bo.le'}
#' {the predefined basic logical operators to compare two elements: ==, >, >=, <, <=.}
#' \item{`fm.buo.neg'}{the predefined basic unary operator for negate.}
#' \item{`fm.buo.sqrt'}{the predefined basic unary operator for square root.}
#' \item{`fm.buo.abs'}{the predefined basic unary operator for absolute value.}
#' \item{`fm.buo.not'}{the predefined logical NOT operator.}
#' \item{`fm.buo.ceil'}{the predefined basic unary operator of computing a ceiling
#' of a number.}
#' \item{`fm.buo.floor'}{the predefined basic unary operator of computing a floor
#' of a number.}
#' \item{`fm.buo.log', `fm.buo.log2' and `fm.buo.log10'}{the predefined basic unary
#' operators of computing log with different bases.}
#' \item{`fm.buo.round'}{the predefined basic unary operator of rounding a value.}
#' \item{`fm.buo.as.int'}{the predefined basic unary operator of casting a numeric
#' value to an integer.}
#' \item{`fm.buo.as.numeric'}{the predefined basic unary operator of casting
#' an integer to a numeric value.}
#' }
#'
#' @param name the name of the basic operator.
#' @return a reference to the specified basic operator.
#' @name fm.basic.op
#' @author Da Zheng <dzheng5@@jhu.edu>
fm.get.basic.op <- function(name)
{
	stopifnot(!is.null(name))
	op <- .Call("R_FM_get_basic_op", name, PACKAGE="FlashR")
	if (!is.null(op))
		new("fm.bo", info=op$info, name=op$name)
}

#' @name fm.basic.op
fm.get.basic.uop <- function(name)
{
	stopifnot(!is.null(name))
	op <- .Call("R_FM_get_basic_uop", name, PACKAGE="FlashR")
	if (!is.null(op))
		new("fm.bo", info=op$info, name=op$name)
}

#' @name fm.basic.op
fm.init.basic.op <- function()
{
	fm.bo.add <<- fm.get.basic.op("add")
	stopifnot(!is.null(fm.bo.add))
	fm.bo.sub <<- fm.get.basic.op("sub")
	stopifnot(!is.null(fm.bo.sub))
	fm.bo.mul <<- fm.get.basic.op("mul")
	stopifnot(!is.null(fm.bo.mul))
	fm.bo.div <<- fm.get.basic.op("div")
	stopifnot(!is.null(fm.bo.div))
	fm.bo.min <<- fm.get.basic.op("min")
	stopifnot(!is.null(fm.bo.min))
	fm.bo.max <<- fm.get.basic.op("max")
	stopifnot(!is.null(fm.bo.max))
	fm.bo.pow <<- fm.get.basic.op("pow")
	stopifnot(!is.null(fm.bo.pow))
	fm.bo.eq <<- fm.get.basic.op("eq")
	stopifnot(!is.null(fm.bo.eq))
	fm.bo.neq <<- fm.get.basic.op("neq")
	stopifnot(!is.null(fm.bo.neq))
	fm.bo.gt <<- fm.get.basic.op("gt")
	stopifnot(!is.null(fm.bo.gt))
	fm.bo.ge <<- fm.get.basic.op("ge")
	stopifnot(!is.null(fm.bo.ge))
	fm.bo.lt <<- fm.get.basic.op("lt")
	stopifnot(!is.null(fm.bo.lt))
	fm.bo.le <<- fm.get.basic.op("le")
	stopifnot(!is.null(fm.bo.le))
	fm.bo.or <<- fm.get.basic.op("|")
	stopifnot(!is.null(fm.bo.or))
	fm.bo.and <<- fm.get.basic.op("&")
	stopifnot(!is.null(fm.bo.and))

	fm.bo.count <<- fm.get.basic.op("count")
	stopifnot(!is.null(fm.bo.count))
	fm.bo.which.max <<- fm.get.basic.op("which.max")
	stopifnot(!is.null(fm.bo.which.max))
	fm.bo.which.min <<- fm.get.basic.op("which.min")
	stopifnot(!is.null(fm.bo.which.min))
	fm.bo.euclidean <<- fm.get.basic.op("euclidean")
	stopifnot(!is.null(fm.bo.euclidean))

	fm.buo.neg <<- fm.get.basic.uop("neg")
	stopifnot(!is.null(fm.buo.neg))
	fm.buo.sqrt <<- fm.get.basic.uop("sqrt")
	stopifnot(!is.null(fm.buo.sqrt))
	fm.buo.abs <<- fm.get.basic.uop("abs")
	stopifnot(!is.null(fm.buo.abs))
	fm.buo.not <<- fm.get.basic.uop("not")
	stopifnot(!is.null(fm.buo.not))
	fm.buo.ceil <<- fm.get.basic.uop("ceil")
	stopifnot(!is.null(fm.buo.ceil))
	fm.buo.floor <<- fm.get.basic.uop("floor")
	stopifnot(!is.null(fm.buo.floor))
	fm.buo.log <<- fm.get.basic.uop("log")
	stopifnot(!is.null(fm.buo.log))
	fm.buo.log2 <<- fm.get.basic.uop("log2")
	stopifnot(!is.null(fm.buo.log2))
	fm.buo.log10 <<- fm.get.basic.uop("log10")
	stopifnot(!is.null(fm.buo.log10))
	fm.buo.round <<- fm.get.basic.uop("round")
	stopifnot(!is.null(fm.buo.round))
	fm.buo.as.int <<- fm.get.basic.uop("as.int")
	stopifnot(!is.null(fm.buo.as.int))
	fm.buo.as.numeric <<- fm.get.basic.uop("as.numeric")
	stopifnot(!is.null(fm.buo.as.numeric))
}

fm.create.agg.op <- function(agg, combine, name)
{
	stopifnot(class(agg) == "fm.bo")
	# It's OK if combine doesn't exist.
	if (is.null(combine)) {
		invalid <- c(as.integer(-1), as.integer(0))
		new("fm.agg.op", agg=agg@info, combine=invalid, name=name)
	}
	else {
		stopifnot(class(combine) == "fm.bo")
		new("fm.agg.op", agg=agg@info, combine=combine@info, name=name)
	}
}

#' Aggregation on a FlashMatrix object.
#'
#' This function accepts a basic operator and perform aggregation on
#' the FlashMatrix object with the basic operator.
#'
#' `fm.agg' aggregates over the entire object.
#'
#' `fm.agg.lazy' aggregates over the entire object, but it performs aggregation
#' lazily.
#'
#' `fm.agg.mat' aggregates on the rows or columns of a matrix. It performs
#' aggregation on the shorter dimension lazily, but on the longer dimension
#' immediately.
#'
#' `fm.agg.mat.lazy' aggregates on the rows or columns of a matrix and performs
#' aggregation lazily regardless the dimension.
#'
#' @param fm a FlashMatrix object
#' @param op the reference or the name of a predefined basic operator or
#' the reference to an aggregation operator returned by `fm.create.agg.op'.
#' @param margin the subscript which the function will be applied over.
#' @return `fm.agg' returns a scalar, `fm.agg.mat' returns a FlashMatrix vector,
#' `fm.agg.lazy' and `fm.agg.mat.lazy' return a FlashMatrix sink matrix.
#' @name fm.basic.op
fm.agg <- function(fm, op, test.na)
{
	stopifnot(!is.null(fm) && !is.null(op))
	stopifnot(class(fm) == "fmV" || class(fm) == "fm")
	if (class(op) == "character")
		op <- fm.get.basic.op(op)
	if (class(op) == "fm.bo")
		op <- fm.create.agg.op(op, op, op@name)
	stopifnot(class(op) == "fm.agg.op")
	.Call("R_FM_agg", fm, op, PACKAGE="FlashR")
}

#' @name fm.basic.op
fm.agg.lazy <- function(fm, op)
{
	stopifnot(!is.null(fm) && !is.null(op))
	stopifnot(class(fm) == "fmV" || class(fm) == "fm")
	if (class(op) == "character")
		op <- fm.get.basic.op(op)
	if (class(op) == "fm.bo")
		op <- fm.create.agg.op(op, op, op@name)
	stopifnot(class(op) == "fm.agg.op")
	ret <- .Call("R_FM_agg_lazy", fm, op, PACKAGE="FlashR")
	new.fmSink(ret)
}

#' @name fm.basic.op
fm.agg.mat <- function(fm, margin, op, test.na)
{
	stopifnot(!is.null(fm) && !is.null(op))
	stopifnot(class(fm) == "fm")
	if (class(op) == "character")
		op <- fm.get.basic.op(op)
	if (class(op) == "fm.bo")
		op <- fm.create.agg.op(op, op, op@name)
	stopifnot(class(op) == "fm.agg.op")
	ret <- .Call("R_FM_agg_mat", fm, as.integer(margin), op, PACKAGE="FlashR")
	new.fmV(ret)
}

#' @name fm.basic.op
fm.agg.mat.lazy <- function(fm, margin, op)
{
	stopifnot(!is.null(fm) && !is.null(op))
	stopifnot(class(fm) == "fm")
	if (class(op) == "character")
		op <- fm.get.basic.op(op)
	if (class(op) == "fm.bo")
		op <- fm.create.agg.op(op, op, op@name)
	stopifnot(class(op) == "fm.agg.op")
	ret <- .Call("R_FM_agg_mat_lazy", fm, as.integer(margin), op,
				 PACKAGE="FlashR")
	new.fmSink(ret)
}

fm.env <- new.env()
fm.env$fm.test.na <- TRUE

fm.set.test.na <- function(val)
{
	fm.env$fm.test.na <- val
}

fm.set.na <- function(in1, in2, res)
{
	# This is a special function that helps to test and set NA on
	# the computation results, so it shouldn't call any functions that
	# try to test and set NA on the result.
	if (is.null(res))
		return(NULL)
	else if (!fm.env$fm.test.na)
		return(res)
	else if (fm.typeof(res) == "logical")
		# is.na always return TRUE or FALSE, we don't need to test and set NA
		# on the result. If we do, we'll get infinite recursive calls.
		ifelse(fm.mapply2(is.na(in1), is.na(in2), fm.bo.or,
						  FALSE), NA, res)
	else if (fm.typeof(res) == "integer")
		ifelse(fm.mapply2(is.na(in1), is.na(in2), fm.bo.or,
						  FALSE), as.integer(NA), res)
	else if (fm.typeof(res) == "double")
		ifelse(fm.mapply2(fm.is.na.only(in1), fm.is.na.only(in2),
						  fm.bo.or, FALSE), as.double(NA), res)
	else
		# In this case, we don't do anything with the result.
		res
}

fm.mapply2.fm <- function(o1, o2, FUN, set.na=TRUE)
{
	if (class(FUN) == "character")
		FUN <- fm.get.basic.op(FUN)
	stopifnot(class(FUN) == "fm.bo")
	stopifnot(dim(o1)[2] == dim(o2)[2] && dim(o1)[1] == dim(o2)[1])
	ret <- .Call("R_FM_mapply2", FUN, o1, o2, PACKAGE="FlashR")
	ret <- new.fm(ret)
	if (set.na)
		ret <- fm.set.na(o1, o2, ret)
	ret
}

fm.mapply2.fmV <- function(o1, o2, FUN, set.na=TRUE)
{
	if (class(FUN) == "character")
		FUN <- fm.get.basic.op(FUN)
	stopifnot(class(FUN) == "fm.bo")
	stopifnot(length(o1) == length(o2))
	ret <- .Call("R_FM_mapply2", FUN, o1, o2, PACKAGE="FlashR")
	ret <- new.fmV(ret)
	if (set.na)
		ret <- fm.set.na(o1, o2, ret)
	ret
}

fm.mapply2.fm.m <- function(o1, o2, FUN, set.na=TRUE)
{
	if (class(FUN) == "character")
		FUN <- fm.get.basic.op(FUN)
	stopifnot(class(FUN) == "fm.bo")
	stopifnot(nrow(o1) == nrow(o2))
	stopifnot(ncol(o1) == ncol(o2))
	o2 <- fm.conv.R2FM(o2)
	if (is.null(o2))
		return(NULL)
	ret <- .Call("R_FM_mapply2", FUN, o1, o2, PACKAGE="FlashR")
	ret <- new.fm(ret)
	if (set.na)
		ret <- fm.set.na(o1, o2, ret)
	ret
}

fm.mapply2.m.fm <- function(o1, o2, FUN, set.na=TRUE)
{
	if (class(FUN) == "character")
		FUN <- fm.get.basic.op(FUN)
	stopifnot(class(FUN) == "fm.bo")
	stopifnot(nrow(o1) == nrow(o2))
	stopifnot(ncol(o1) == ncol(o2))
	o1 <- fm.conv.R2FM(o1)
	if (is.null(o1))
		return(NULL)
	ret <- .Call("R_FM_mapply2", FUN, o1, o2, PACKAGE="FlashR")
	ret <- new.fm(ret)
	if (set.na)
		ret <- fm.set.na(o1, o2, ret)
	ret
}

fm.mapply2.fm.ANY <- function(o1, o2, FUN, set.na=TRUE)
{
	if (class(FUN) == "character")
		FUN <- fm.get.basic.op(FUN)
	stopifnot(class(FUN) == "fm.bo")
	stopifnot(is.vector(o2) || fm.is.vector(o2))
	if (length(o2) > 1) {
		o2 <- fm.conv.R2FM(o2)
		if (is.null(o2))
			return(NULL)
		fm.mapply.col(o1, o2, FUN, set.na)
	}
	else {
		ret <- .Call("R_FM_mapply2_AE", FUN, o1, o2, PACKAGE="FlashR")
		ret <- new.fm(ret)
		if (set.na)
			ret <- fm.set.na(o1, o2, ret)
		ret
	}
}

fm.mapply2.ANY.fm <- function(o1, o2, FUN, set.na=TRUE)
{
	if (class(FUN) == "character")
		FUN <- fm.get.basic.op(FUN)
	stopifnot(class(FUN) == "fm.bo")
	stopifnot(is.vector(o1) || fm.is.vector(o1))
	if (length(o1) > 1) {
		print("don't support this operation yet.")
		NULL
	}
	else {
		ret <- .Call("R_FM_mapply2_EA", FUN, o1, o2, PACKAGE="FlashR")
		ret <- new.fm(ret)
		if (set.na)
			ret <- fm.set.na(o1, o2, ret)
		ret
	}
}

fm.mapply2.fmV.ANY <- function(o1, o2, FUN, set.na=TRUE)
{
	if (class(FUN) == "character")
		FUN <- fm.get.basic.op(FUN)
	stopifnot(class(FUN) == "fm.bo")
	stopifnot(is.vector(o2) || fm.is.vector(o2))
	if (length(o2) > 1) {
		stopifnot(length(o2) == length(o1))
		o2 <- fm.conv.R2FM(o2)
		if (is.null(o2))
			return(NULL)
		ret <- .Call("R_FM_mapply2", FUN, o1, o2, PACKAGE="FlashR")
	}
	else
		ret <- .Call("R_FM_mapply2_AE", FUN, o1, o2, PACKAGE="FlashR")
	ret <- new.fmV(ret)
	if (set.na)
		ret <- fm.set.na(o1, o2, ret)
	ret
}

fm.mapply2.ANY.fmV <- function(o1, o2, FUN, set.na=TRUE)
{
	if (class(FUN) == "character")
		FUN <- fm.get.basic.op(FUN)
	stopifnot(class(FUN) == "fm.bo")
	stopifnot(is.vector(o1) || fm.is.vector(o1))
	if (length(o1) == 1)
		ret <- .Call("R_FM_mapply2_EA", FUN, o1, o2, PACKAGE="FlashR")
	else {
		stopifnot(length(o1) == length(o2))
		o1 <- fm.conv.R2FM(o1)
		if (is.null(o1))
			return(NULL)
		ret <- .Call("R_FM_mapply2", FUN, o1, o2, PACKAGE="FlashR")
	}
	ret <- new.fmV(ret)
	if (set.na)
		ret <- fm.set.na(o1, o2, ret)
	ret
}

fm.mapply.row <- function(o1, o2, FUN, set.na=TRUE)
{
	stopifnot(class(o1) == "fm" && class(o2) == "fmV")
	if (class(FUN) == "character")
		FUN <- fm.get.basic.op(FUN)
	stopifnot(class(FUN) == "fm.bo")
	ret <- .Call("R_FM_mapply2_MV", o1, o2, as.integer(1), FUN,
				 PACKAGE="FlashR")
	ret <- new.fm(ret)
	if (set.na)
		ret <- fm.set.na(o1, o2, ret)
	ret
}

fm.mapply.col <- function(o1, o2, FUN, set.na=TRUE)
{
	stopifnot(class(o1) == "fm" && class(o2) == "fmV")
	if (class(FUN) == "character")
		FUN <- fm.get.basic.op(FUN)
	stopifnot(class(FUN) == "fm.bo")
	ret <- .Call("R_FM_mapply2_MV", o1, o2, as.integer(2), FUN,
				 PACKAGE="FlashR")
	ret <- new.fm(ret)
	if (set.na)
		ret <- fm.set.na(o1, o2, ret)
	ret
}

#' Apply a Function to two FlashMatrix vectors/matrices.
#'
#' `mapply' applies `FUN' to the first elements of each vector, the second
#' elements, the third elements, and so on. Two vectors/matrices should have
#' the same shape. Currently, `mapply' only accepts predefined basic operators
#' returned by `fm.get.basic.op'.
#'
#' @param o1, o2 a FlashMatrix vector/matrix.
#' @param FUN the reference or the name of one of the predefined basic binary
#' operators.
#' @return a FlashMatrix vector/matrix.
#' @name fm.mapply2
#' @author Da Zheng <dzheng5@@jhu.edu>
setGeneric("fm.mapply2", function(o1, o2, FUN, set.na) standardGeneric("fm.mapply2"))
setMethod("fm.mapply2",
		  signature(o1 = "fm", o2 = "fm", FUN = "ANY", set.na="logical"),
		  fm.mapply2.fm)
setMethod("fm.mapply2",
		  signature(o1 = "fmV", o2 = "fmV", FUN = "ANY", set.na="logical"),
		  fm.mapply2.fmV)
setMethod("fm.mapply2",
		  signature(o1 = "fm", o2 = "fmV", FUN = "ANY", set.na="logical"),
		  fm.mapply.col)
setMethod("fm.mapply2",
		  signature(o1 = "fmV", o2 = "fm", FUN = "ANY", set.na="logical"),
		  function(o1, o2, FUN, set.na) {
			  print("This isn't supported currently.")
			  NULL
		  })
setMethod("fm.mapply2",
		  signature(o1 = "fm", o2 = "matrix", FUN = "ANY", set.na="logical"),
		  fm.mapply2.fm.m)
setMethod("fm.mapply2",
		  signature(o1 = "matrix", o2 = "fm", FUN = "ANY", set.na="logical"),
		  fm.mapply2.m.fm)
setMethod("fm.mapply2",
		  signature(o1 = "fm", o2 = "ANY", FUN = "ANY", set.na="logical"),
		  fm.mapply2.fm.ANY)
setMethod("fm.mapply2",
		  signature(o1 = "ANY", o2 = "fm", FUN = "ANY", set.na="logical"),
		  fm.mapply2.ANY.fm)
setMethod("fm.mapply2",
		  signature(o1 = "fmV", o2 = "ANY", FUN = "ANY", set.na="logical"),
		  fm.mapply2.fmV.ANY)
setMethod("fm.mapply2",
		  signature(o1 = "ANY", o2 = "fmV", FUN = "ANY", set.na="logical"),
		  fm.mapply2.ANY.fmV)

fm.set.na1 <- function(input, res)
{
	# This is a special function that helps to test and set NA on
	# the computation results, so it shouldn't call any functions that
	# try to test and set NA on the result.
	if (is.null(res))
		return(NULL)
	else if (!fm.env$fm.test.na)
		return(res)
	else if (fm.typeof(res) == "logical")
		# is.na always return TRUE or FALSE, we don't need to test and set NA
		# on the result. If we do, we'll get infinite recursive calls.
		ifelse(is.na(input), NA, res)
	else if (fm.typeof(res) == "integer")
		ifelse(is.na(input), as.integer(NA), res)
	else if (fm.typeof(res) == "double") {
		ifelse(fm.is.na.only(input), as.double(NA), res)
	}
	else
		# In this case, we don't do anything with the result.
		res
}

fm.sapply.fm <- function(o, FUN, set.na=TRUE)
{
	if (class(FUN) == "character")
		FUN <- fm.get.basic.uop(FUN)
	ret <- .Call("R_FM_sapply", FUN, o, PACKAGE="FlashR")
	ret <- new.fm(ret)
	if (set.na)
		ret <- fm.set.na1(o, ret)
	ret
}

fm.sapply.fmV <- function(o, FUN, set.na=TRUE)
{
	if (class(FUN) == "character")
		FUN <- fm.get.basic.uop(FUN)
	ret <- .Call("R_FM_sapply", FUN, o, PACKAGE="FlashR")
	ret <- new.fmV(ret)
	if (set.na)
		ret <- fm.set.na1(o, ret)
	ret
}

#' Apply a Function to a FlashMatrix vector/matrix.
#'
#' `sapply' applies `FUN' to every element of a vector/matrix.
#' Currently, `sapply' only accepts predefined basic operators
#' returned by `fm.get.basic.uop'.
#'
#' @param o a FlashMatrix vector/matrix.
#' @param FUN the reference or the name of a predefined uniary operator.
#' @return a FlashMatrix vector/matrix.
#' @name fm.sapply
#' @author Da Zheng <dzheng5@@jhu.edu>
setGeneric("fm.sapply", function(o, FUN, set.na)  standardGeneric("fm.sapply"))
setMethod("fm.sapply", signature(o = "fm", FUN = "ANY", set.na="logical"),
		  fm.sapply.fm)
setMethod("fm.sapply", signature(o = "fmV", FUN = "ANY", set.na="logical"),
		  fm.sapply.fmV)

#' Groupby on a FlashMatrix vector.
#'
#' `fm.sgroupby' groups elements in a vector based on corresponding `labels'
#' and applies `FUN' to the elements in each group. `FUN' is an aggregation
#' operator.
#'
#' `fm.groupby' groups rows/columns of a matrix based on corresponding `labels'
#' and applies `FUN' to the rows/columns in each group. `FUN' is an aggregation
#' operator.
#'
#' @param obj a FlashMatrix vector or matrix
#' @param margin the subscript which the function will be applied over.
#' E.g., for a matrix, `1' indicates rows, `2' indicates columns.
#' @param labels a FlashMatrix vector that indicates the labels of
#' each element in `obj'.
#' @param FUN an aggregation operator returned by `fm.create.agg.op'.
#' @return `fm.sgroupby' returns a data frame, where the column `val' stores
#' all of the unique values in the original data container, and the column
#' `agg' stores the aggregate result of the corresponding value.
#' @name fm.groupby
#' @author Da Zheng <dzheng5@@jhu.edu>
fm.sgroupby <- function(obj, FUN)
{
	stopifnot(class(obj) == "fmV" || class(obj) == "fmFactorV")
	if (class(FUN) == "character")
		FUN <- fm.get.basic.op(FUN)
	if (class(FUN) == "fm.bo")
		FUN <- fm.create.agg.op(FUN, FUN, FUN@name)
	stopifnot(class(FUN) == "fm.agg.op")
	res <- .Call("R_FM_sgroupby", obj, FUN, PACKAGE="FlashR")
	if (is.null(res))
		return(NULL)
	list(val=new.fmV(res$val), Freq=new.fmV(res$agg))
}

#' @name fm.groupby
fm.groupby <- function(fm, margin, factor, FUN)
{
	stopifnot(class(fm) == "fm")
	if (class(FUN) == "character")
		FUN <- fm.get.basic.op(FUN)
	if (class(FUN) == "fm.bo")
		FUN <- fm.create.agg.op(FUN, FUN, FUN@name)
	stopifnot(class(FUN) == "fm.agg.op")
	stopifnot(class(factor) == "fmFactorV")
	orig.margin <- margin
	if (margin == 1) {
		margin <- 2
		fm <- t(fm)
	}
	res <- .Call("R_FM_groupby", fm, as.integer(margin), factor, FUN,
				 PACKAGE="FlashR")
	res <- new.fm(res)
	if (orig.margin == 1)
		t(res)
	else
		res
}

#' Transpose a FlashMatrix matrix.
#'
#' @param m a FlashMatrix matrix
#' @return a FlashMatrix matrix
#' @author Da Zheng <dzheng5@@jhu.edu>
fm.t <- function(m)
{
	stopifnot(!is.null(m))
	stopifnot(class(m) == "fm")
	ret <- .Call("R_FM_transpose", m, PACKAGE="FlashR")
	new.fm(ret)
}

#' Set the specified column of a FlashMatrix matrix.
#'
#' @param fm A flashMatrixR matrix
#' @param idxs an array of column indices in fm.
#' @param m2 A flashMatrixR matrix or vector.
#' @return a logical value that indicates whether the operation is successful.
#' @author Da Zheng <dzheng5@@jhu.edu>
fm.set.cols <- function(fm, idxs, m2)
{
	print("fm.set.cols isn't supported right now.")
	stopifnot(FALSE)
#	stopifnot(!is.null(fm) && !is.null(idxs) && !is.null(m2))
#	stopifnot(class(fm) == "fm")
#	stopifnot(class(m2) == "fmV" || class(m2) == "fm")
#	.Call("R_FM_set_cols", fm, as.integer(idxs), m2, PACKAGE="FlashR")
}

#' Get a submatrix from a FlashMatrix matrix
#'
#' `fm.get.rows' gets specified rows in a FM matrix.
#' `fm.get.cols' gets specified columns in a FM matrix.
#'
#' @param fm A FlashMatrix matrix
#' @param idxs an array of column indices in fm.
#' @return a submatrix
#' @name fm.get.eles
#' @author Da Zheng <dzheng5@@jhu.edu>
fm.get.cols <- function(fm, idxs)
{
	stopifnot(!is.null(fm) && !is.null(idxs))
	stopifnot(class(fm) == "fm")
	ret <- .Call("R_FM_get_submat", fm, as.integer(2), as.integer(idxs),
				 PACKAGE="FlashR")
	if (!is.null(ret))
		new.fm(ret)
	else
		NULL
}

#' @rdname fm.get.eles
fm.get.rows <- function(fm, idxs)
{
	stopifnot(!is.null(fm) && !is.null(idxs))
	stopifnot(class(fm) == "fm")
	ret <- .Call("R_FM_get_submat", fm, as.integer(1), as.integer(idxs),
				 PACKAGE="FlashR")
	if (!is.null(ret))
		new.fm(ret)
	else
		NULL
}

fm.set.materialize.level <- function(fm, level)
{
	stopifnot(!is.null(fm))
	stopifnot(class(fm) == "fm" || class(fm) == "fmV")
	.Call("R_FM_set_materialize_level", fm, as.integer(level), PACKAGE="FlashR")
}

fm.materialize <- function(...)
{
	args <- list(...)
	if (length(args) == 0)
		stop("no arguments")
	else if (length(args) == 1) {
		obj <- args[[1]]
		stopifnot(!is.null(obj))
		stopifnot(class(obj) == "fm" || class(obj) == "fmV"
				  || class(obj) == "fmSink")
		ret <- .Call("R_FM_materialize", obj, PACKAGE="FlashR")
		if (class(obj) == "fm")
			new.fm(ret)
		else
			new.fmV(ret)
	}
	else {
		for (obj in args) {
			stopifnot(!is.null(obj))
			stopifnot(class(obj) == "fm" || class(obj) == "fmV"
					  || class(obj) == "fmSink")
		}
		rets <- .Call("R_FM_materialize_list", args, PACKAGE="FlashR")
		stopifnot(length(rets) == length(args))
		for (i in 1:length(args)) {
			if (class(args[[i]]) == "fm")
				rets[[i]] <- new.fm(rets[[i]])
			else
				rets[[i]] <- new.fmV(rets[[i]])
		}
		rets
	}
}

#' Write a FlashMatrix object (vector/matrix) to a file
#'
#' @param fm a FlashMatrix object.
#' @param file a file in the local filesystem.
#' @return a logical value. True if the object is written to a file
#' successfully. Otherwise, FALSE.
#' @author Da Zheng <dzheng5@@jhu.edu>
fm.write.obj <- function(fm, file)
{
	stopifnot(!is.null(fm))
	stopifnot(class(fm) == "fm" || class(fm) == "fmV")
	.Call("R_FM_write_obj", fm, file, PACKAGE="FlashR")
}

#' Read a FlashMatrix object (vector/matrix) from a file.
#'
#' @param file a file in the local filesystem.
#' @return a FlashMatrix object (vector/matrix)
#' @author Da Zheng <dzheng5@@jhu.edu>
fm.read.obj <- function(file)
{
	ret <- .Call("R_FM_read_obj", file, PACKAGE="FlashR")
	if (class(ret) == "fmV")
		new.fmV(ret)
	else
		new.fm(ret)
}

#' Eigensolver
#'
#' Compute eigenvalues/vectors of the adjacency matrix of an undirected graph.
#'
#' This eigensolver is powered by Anasazi package of Trilinos.
#'
#' @param func The function to perform the matrix-vector multiplication.
#' @param extra Extra argument to supply to `func'.
#' @param sym Logical scalar, whether the input matrix is symmetric.
#' @param options Options to Anasazi.
#'
#' The `options' argument specifies what kind of computation to perform.
#' It is a list with the following members, which correspond directly to
#' Anasazi parameters:
#'
#' nev Numeric scalar. The number of eigenvalues to be computed.
#'
#' solver String. The name of the eigensolver to solve the eigenproblems.
#'				Currently, it supports three eigensolvers: KrylovSchur,
#'				Davidson and LOBPCG. KrylovSchur is the default eigensolver.
#'
#' tol Numeric scalar. Stopping criterion: the relative accuracy of
#'				the Ritz value is considered acceptable if its error is less
#'				than `tol' times its estimated value.
#'
#' block_size Numeric scalar. The eigensolvers use a block extension of an
#'				eigensolver algorithm. The block size determines the number
#'				of the vectors that operate together.
#'
#' num_blocks Numeric scalar. The number of blocks to compute eigenpairs.
#'
#' which Specify which eigenvalues/vectors to compute, character
#'              constant with exactly two characters.
#' Possible values for symmetric input matrices:
#' \itemize{
#' \item{"LA"}{Compute `nev' largest (algebraic) eigenvalues.}
#' \item{"SA"}{Compute `nev' smallest (algebraic) eigenvalues.}
#' \item{"LM"}{Compute `nev' largest (in magnitude) eigenvalues.}
#' \item{"SM"}{Compute `nev' smallest (in magnitude) eigenvalues.}
#' }
#'
#' @return A named list with the following members:
#'         vals: Numeric vector, the desired eigenvalues.
#'         vecs: Numeric matrix, the desired eigenvectors as columns.
#' @name fm.eigen
#' @author Da Zheng <dzheng5@@jhu.edu>
fm.eigen <- function(func, extra=NULL, sym=TRUE, options=NULL,
					 env = parent.frame())
{
	if (!sym) {
		print("fm.eigen only supports symmetric matrices")
		return(NULL)
	}
	if (is.loaded("R_FM_eigen")) {
		ret <- .Call("R_FM_eigen", as.function(func), extra, as.logical(sym),
					 as.list(options), PACKAGE="FlashR")
		ret$vecs <- new.fm(ret$vecs)
		ret
	}
	else {
		arpack.opts <- arpack_defaults
		if (!is.null(options)) {
			if (!is.null(options[["n"]]))
				arpack.opts$n <- options$n
			if (!is.null(options[["nev"]]))
				arpack.opts$nev <- options$nev
			if (!is.null(options[["tol"]]))
				arpack.opts$tol <- options$tol
			if (!is.null(options[["block_size"]])
				&& !is.null(options[["num_blocks"]]))
				arpack.opts$ncv <- options$block_size * options$num_blocks
			# If we compute a small number of eigenvalues, we need a larger
			# subspace.
			else if (arpack.opts$nev <= 2)
				arpack.opts$ncv <- arpack.opts$nev * 2 + 3
			else
				arpack.opts$ncv <- arpack.opts$nev * 2
			if (!is.null(options[["which"]]))
				arpack.opts$which <- options$which
		}
		arpack.fun <- function(x, extra) {
			# TODO this might be slow
			fm.x <- fm.as.vector(x)
			ret <- func(fm.x, extra)
			# TODO this might be slow
			as.vector(ret)
		}
		arpack.ret <- arpack(arpack.fun, extra, sym, arpack.opts, env)
		list(vals=arpack.ret$values, vecs=fm.as.matrix(arpack.ret$vectors),
			 options=arpack.ret$options)
	}
}

#' Combine FlashR matrices by rows or columns.
#'
#' Take a list of FlashR matrices and combine them by Columns or rows
#' respectively.
#'
#' @param ... A list of FlashR matrices.
#'
#' The function of these two functions is similar to the R counterparts:
#' `rbind' and `cbind'. ' Currently, `fm.rbind' only supports combining
#' a list of wide matrices. `fm.cbind' can only combine a list of tall matrices.
#'
#' @return A FlashR matrix.
#' @name fm.bind
#' @author Da Zheng <dzheng5@@jhu.edu>
fm.rbind <- function(...)
{
	args <- list(...)
	nargs <- length(args)
	if (nargs < 2) {
		return(args[1])
	}
	for (fm in args) {
		if (class(fm) != "fm") {
			print("fm.rbind only works on FlashMatrix matrix")
			return(NULL)
		}
	}
	ret <- .Call("R_FM_rbind", args, PACKAGE="FlashR")
	new.fm(ret)
}

#' @name fm.bind
fm.cbind <- function(...)
{
	args <- list(...)
	nargs <- length(args)
	if (nargs < 2) {
		return(args[1])
	}
	for (fm in args) {
		if (class(fm) != "fm") {
			print("fm.cbind only works on FlashMatrix matrix")
			return(NULL)
		}
	}
	ret <- .Call("R_FM_cbind", args, PACKAGE="FlashR")
	new.fm(ret)
}

setMethod("ifelse", signature(test = "fm", yes = "fm", no = "ANY"),
		  function(test, yes, no) {
			  ret <- .Call("R_FM_ifelse_no", test, yes, no, PACKAGE="FlashR")
			  new.fm(ret)
		  })

setMethod("ifelse", signature(test = "fm", yes = "ANY", no = "fm"),
		  function(test, yes, no) {
			  ret <- .Call("R_FM_ifelse_yes", test, yes, no, PACKAGE="FlashR")
			  new.fm(ret)
		  })

setMethod("ifelse", signature(test = "fmV", yes = "fmV", no = "ANY"),
		  function(test, yes, no) {
			  ret <- .Call("R_FM_ifelse_no", test, yes, no, PACKAGE="FlashR")
			  new.fmV(ret)
		  })

setMethod("ifelse", signature(test = "fmV", yes = "ANY", no = "fmV"),
		  function(test, yes, no) {
			  ret <- .Call("R_FM_ifelse_yes", test, yes, no, PACKAGE="FlashR")
			  new.fmV(ret)
		  })

# This returns NA of the right type.
get.na <- function(type) {
	if (type == "logical") {
		NA
	}
	else if (type == "integer") {
		as.integer(NA)
	}
	else {
		stop("unsupported type for NA")
	}
}

fm.is.na.only <- function(fm)
{
	stopifnot(class(fm) == "fm" || class(fm) == "fmV")
	ret <- .Call("R_FM_isna", fm, TRUE, PACKAGE="FlashR")
	if (class(fm) == "fm")
		new.fm(ret)
	else
		new.fmV(ret)
}

# These functions shouldn't call any functions that try to test and
# set NA on the result.

# is.na for float-point needs to be handled differently because C/C++
# interprets NA as NaN.
setMethod("is.na", signature(x = "fm"), function(x) {
		  if (typeof(x) == "double") {
			  ret <- .Call("R_FM_isna", x, FALSE, PACKAGE="FlashR")
			  new.fm(ret)
		  }
		  else
			  # The result has to be TRUE or FALSE.
			  fm.mapply2(x, get.na(typeof(x)), fm.bo.eq, FALSE)
		  })
setMethod("is.na", signature(x = "fmV"), function(x) {
		  if (typeof(x) == "double") {
			  ret <- .Call("R_FM_isna", x, FALSE, PACKAGE="FlashR")
			  new.fmV(ret)
		  }
		  else
			  # The result has to be TRUE or FALSE.
			  fm.mapply2(x, get.na(typeof(x)), fm.bo.eq, FALSE)
		  })

setMethod("is.nan", signature(x = "fm"), function(x) {
		  if (typeof(x) == "double") {
			  ret <- .Call("R_FM_isnan", x, PACKAGE="FlashR")
			  new.fm(ret)
		  }
		  else
			  # TODO construct this matrix shouldn't have data movement.
			  fm.matrix(fm.rep.int(FALSE, length(x)), nrow(x), ncol(x))
		  })
setMethod("is.nan", signature(x = "fmV"), function(x) {
		  if (typeof(x) == "double") {
			  ret <- .Call("R_FM_isnan", x, PACKAGE="FlashR")
			  new.fmV(ret)
		  }
		  else
			  fm.rep.int(FALSE, length(x))
		  })

setMethod("is.infinite", signature(x = "fm"), function(x) {
		  if (typeof(x) == "double")
			  fm.mapply2(x, Inf, fm.bo.eq, FALSE)
		  else
			  # TODO construct this matrix shouldn't have data movement.
			  fm.matrix(fm.rep.int(FALSE, length(x)), nrow(x), ncol(x))
		  })
setMethod("is.infinite", signature(x = "fmV"), function(x) {
		  if (typeof(x) == "double")
			  fm.mapply2(x, Inf, fm.bo.eq, FALSE)
		  else
			  fm.rep.int(FALSE, length(x))
		  })

fm.test.na.finite <- function(input, res)
{
	if (is.null(res))
		return(NULL)
	else {
		# The result must be float points.
		ifelse(fm.mapply2(is.na(input), is.nan(input), fm.bo.or,
						  FALSE), FALSE, res)
	}
}
setMethod("is.finite", signature(x = "fm"), function(x) {
		  if (typeof(x) == "double")
			  fm.test.na.finite(x, fm.mapply2(x, Inf, fm.bo.neq, FALSE))
		  else
			  # TODO construct this matrix shouldn't have data movement.
			  fm.matrix(fm.rep.int(FALSE, length(x)), nrow(x), ncol(x))
		  })
setMethod("is.finite", signature(x = "fmV"), function(x) {
		  if (typeof(x) == "double")
			  fm.test.na.finite(x, fm.mapply2(x, Inf, fm.bo.neq, FALSE))
		  else
			  fm.rep.int(FALSE, length(x))
		  })

fm.print.mat.info <- function(fm)
{
	stopifnot(class(fm) == "fm" || class(fm) == "fmV")
	ret <- .Call("R_FM_print_mat_info", fm, PACKAGE="FlashR")
}
