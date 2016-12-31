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

#' An S4 class to represent a FlashR matrix.
#'
#' @slot pointer points a matrix object in C.
#' @slot name a string indicating the name of the matrix.
#' @slot nrow a numeric value indicating the number of rows.
#' @slot ncol a numeric value indicating the number of columns.
#' @slot type a string indicating the type of the matrix. e.g., sparse or dense.
#' @slot ele_type a string indicating the element type in the matrix.
#' @slot attrs a list that stores the attributes of the matrix.
setClass("fm", representation(pointer = "externalptr", name = "character",
							  nrow = "numeric", ncol = "numeric",
							  type="character", ele_type="character",
							  attrs="list"))

#' An S4 class to represent a FlashR vector.
#'
#' @slot pointer points a vector object in C.
#' @slot name a string indicating the name of the vector.
#' @slot len a numeric value indicating the number of elements.
#' @slot type a string indicating the type of the vector. This field isn't
#'            used for a vector.
#' @slot ele_type a string indicating the element type in the vector.
#' @slot attrs a list that stores the attributes of the matrix.
setClass("fmV", representation(pointer = "externalptr", name = "character",
							   len = "numeric", type="character",
							   ele_type="character", attrs="list"))

#' An S4 class to represent a FlashR factor vector. It inherits from
#' a FlashR vector.
#'
#' @slot pointer points a vector object in C.
#' @slot name a string indicating the name of the vector.
#' @slot len a numeric value indicating the number of elements.
#' @slot type a string indicating the type of the vector. This field isn't
#'            used for a vector.
#' @slot ele_type a string indicating the element type in the vector.
#' @slot num.levels an integer indicating the number of levels.
setClass("fmFactorV", representation(num.levels = "integer"), contains = "fmV")

#' An S4 class to represent a binary operator used in generalized matrix
#' operations.
#'
#' @slot info an integer indicating the binary operator registered
#'            in the system.
#' @slot name a string indicating the name of the operator.
setClass("fm.bo", representation(info = "integer", name = "character"))

#' An S4 class to represent an aggregation operator used in aggregation
#' operations in FlashR.
#'
#' @slot agg an integer indicating the operator for computing partial
#'           aggregation results.
#' @slot combine an integer indicating the operator for combining partial
#'               aggregation results and computing the final result.
#' @slot name a string indicating the name of the aggregation operation.
setClass("fm.agg.op", representation(agg = "integer", combine = "integer",
									 name = "character"))

#' An S4 class to represent an apply operator on rows/columns in a matrix.
#' This operator is used in \code{fm.apply}.
#'
#' @slot info an integer indicating the apply operator registered in the system.
#' @slot name a string indicating the name of the operator.
setClass("fm.apply.op", representation(info = "integer", name = "character"))

.new.fm <- function(fm)
{
	if (!is.null(fm))
		new("fm", pointer=fm$pointer, name=fm$name, nrow=fm$nrow, ncol=fm$ncol,
			type=fm$type, ele_type=fm$ele_type)
	else
		NULL
}

.new.fmV <- function(fm)
{
	if (!is.null(fm))
		new("fmV", pointer=fm$pointer, name=fm$name, len=fm$len, type=fm$type,
			ele_type=fm$ele_type)
	else
		NULL
}

#' Reconfigure FlashR
#'
#' \code{fm.set.conf} reconfigures FlashR with the settings in
#' the configuration file. \code{fm.print.conf} prints the current configurations.
#' \code{fm.set.log.level} sets what levels of messages should be logged.
#' \code{fm.print.features} prints the features compiled into FlashR.
#'
#' The configuration file contains a list of key-value pairs. Each line in
#' the file is a key-value pair in the form of "key_name=value".
#' @param conf.file The configuration file.
#' @param level a string, including "debug", "info", "warning", "error", "fatal".
#' @name fm.set.conf
#' @author Da Zheng <dzheng5@@jhu.edu>
NULL

#' @rdname fm.set.conf
fm.set.conf <- function(conf.file)
{
	ret <- .Call("R_FM_set_conf", as.character(conf.file), PACKAGE="FlashR")
	stopifnot(ret);
}

#' @rdname fm.set.conf
fm.set.log.level <- function(level)
{
	.Call("R_FM_set_log_level", as.character(level), PACKAGE="FlashR")
}

#' @rdname fm.set.conf
fm.print.conf <- function()
{
	fm.set.log.level("info")
	.Call("R_SAFS_print_conf", PACKAGE="FlashR")
	.Call("R_FM_print_conf", PACKAGE="FlashR")
	fm.set.log.level("warning")
}

#' @rdname fm.set.conf
fm.print.features <- function()
{
	.Call("R_FM_print_features", PACKAGE="FlashR")
}

#' Load a matrix to FlashR.
#'
#' There are many different ways of loading a matrix to FlashR.
#' \code{fm.load.dense.matrix} loads a dense matrix in the text format from
#' the Linux filesystem.
#' \code{fm.load.dense.matrix.bin} loads a dense matrix in the binary format
#' from the Linux filesystem.
#' \code{fm.load.sparse.matrix} loads a FlashR sparse matrix from files.
#' The matrix in the file is in the FlashR format.
#' \code{fm.get.dense.matrix} returns a named dense matrix that has already
#' been loaded to FlashR.
#'
#' If a user provides \code{name} and \code{in.mem} is \code{TRUE}, the created
#' vector/matrix will be kept on disks persistently. That is, even if a user
#' exits from R, the vector/matrix will still be kept on disks. A user can
#' access to the dense matrix with \code{fm.get.dense.matrix} the next time
#' when he/she opens FlashR.
#'
#' @param src.file a string that indicates the file in the Linux filesystem
#'        that stores data to be loaded to FlashR.
#' @param spm The file that stores the sparse matrix.
#' @param spm.idx The file that stores the index of the sparse matrix.
#' @param t.spm The file that stores the transpose of the sparse matrix.
#' @param t.spm.idx The file that stores the index of the transpose of the sparse matrix.
#' @param in.mem Determine the loaded matrix is stored in memory or on SAFS.
#' @param ele.type The element type in a matrix.
#' @param delim The delimiter of separating elements in the text format.
#' @param nrow the number of rows in the binary dense matrix.
#' @param ncol the number of columns in the binary dense matrix.
#' @param byrow a logical value indicating if the data in the binary matrix
#'              is stored by rows.
#' @param name a string indicating the name of the dense matrix after being
#'        loaded to FlashR.
#' @return a FlashR matrix.
#' @name fm.get.matrix
#' @author Da Zheng <dzheng5@@jhu.edu>
NULL

#' @rdname fm.get.matrix
fm.get.dense.matrix <- function(name)
{
	stopifnot(!is.null(name))
	stopifnot(class(name) == "character")
	m <- .Call("R_FM_get_dense_matrix", as.character(name), PACKAGE="FlashR")
	.new.fm(m)
}

#' @rdname fm.get.matrix
fm.load.dense.matrix <- function(src.file, in.mem, ele.type="D", delim=",",
								 ncol=.Machine$integer.max, name="")
{
	stopifnot(!is.null(src.file))
	if (is.double(ncol))
		ncol = as.integer(ncol)
	if (is.character(src.file))
		src.file <- c(src.file)
	m <- .Call("R_FM_load_dense_matrix", as.character(src.file), as.logical(in.mem),
			   as.character(ele.type), as.character(delim),
			   ncol=ncol, as.character(name), PACKAGE="FlashR")
	.new.fm(m)
}

#' @rdname fm.get.matrix
fm.load.dense.matrix.bin <- function(src.file, in.mem, nrow, ncol, byrow, ele.type,
									 name="")
{
	m <- .Call("R_FM_load_dense_matrix_bin", as.character(src.file),
			   as.logical(in.mem), as.double(nrow), as.double(ncol),
			   as.logical(byrow), as.character(ele.type), as.character(name),
			   PACKAGE="FlashR")
	.new.fm(m)
}

#' @rdname fm.get.matrix
fm.load.sparse.matrix <- function(spm, spm.idx, t.spm=NULL, t.spm.idx=NULL, in.mem=TRUE)
{
	if (is.null(t.spm) || is.null(t.spm.idx))
		m <- .Call("R_FM_load_matrix_sym", as.character(spm), as.character(spm.idx),
				   as.logical(in.mem), PACKAGE="FlashR")
	else
		m <- .Call("R_FM_load_matrix_asym", as.character(spm), as.character(spm.idx),
				   as.character(t.spm), as.character(t.spm.idx), as.logical(in.mem),
			  PACKAGE="FlashR")
	.new.fm(m)
}

#' Create a FlashR vector with replicated elements.
#'
#' @param x A constant initial value
#' @param times The length of the vector to be generated.
#' @return A FlashR vector
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
	.new.fmV(vec)
}

#' Sequence Generation
#'
#' \code{fm.seq.int} creates a FlashR vector with a sequence of numbers.
#' \code{fm.seq.matrix} creates a FlashR matrix with a sequence of numbers.
#'
#' @param from the starting value of the sequence.
#' @param to the end value of the sequence.
#' @param by number. increment of the sequence.
#' @param nrow the number of rows in the generated matrix.
#' @param ncol the number of columns in the generated matrix.
#' @param byrow logical. If \code{FALSE} (the default) the matrix is filled by
#'	            columns, otherwise the matrix is filled by rows.
#' @return a FlashR vector or matrix filled with a sequence of numbers
#'         that belongs to [from, to]
#' @name fm.seq
#' @author Da Zheng <dzheng5@@jhu.edu>
NULL

#' @rdname fm.seq
fm.seq.int <- function(from, to, by)
{
	if ((to - from) / by <= 0) {
		print("we can't generate a vector of 0 elements")
		return(NULL)
	}
	vec <- .Call("R_FM_create_seq", from, to, by, PACKAGE="FlashR")
	.new.fmV(vec)
}

#' @rdname fm.seq
fm.seq.matrix <- function(from, to, nrow, ncol, byrow = FALSE)
{
	mat <- .Call("R_FM_create_seq_matrix", from, to, nrow, ncol, byrow,
				 PACKAGE="FlashR")
	.new.fm(mat)
}

#' The Uniform Distribution
#'
#' \code{fm.runif} creates a FlashR vector with uniformly random numbers.
#' \code{fm.runif.matrix} creates a FlashR matrix with uniformly random
#' numbers.
#'
#' If a user provides \code{name} and \code{in.mem} is \code{TRUE}, the created
#' vector/matrix will be kept on disks persistently. That is, even if a user
#' exits from R, the vector/matrix will still be kept on disks.
#'
#' @param n the number of random numbers to be generated.
#' @param nrow the number of rows in the generated matrix.
#' @param ncol the number of columns in the generated matrix.
#' @param min lower limits of the distribution.
#' @param max upper limits of the distribution.
#' @param in.mem whether the vector is stored in memory.
#' @param name the name of the vector. It's stored on disks, it's used as
#'             the file name.
#' @name fm.runif
#' @author Da Zheng <dzheng5@@jhu.edu>
NULL

#' @rdname fm.runif
fm.runif <- function(n, min=0, max=1, in.mem=TRUE, name="")
{
	n <- floor(as.numeric(n))
	if (n <= 0) {
		print("we can't generate a vector of 0 elements")
		return(NULL)
	}
	vec <- .Call("R_FM_create_randmat", "uniform", n, 1,
				 as.logical(in.mem), as.character(name),
				 list(min=as.double(min), max=as.double(max)), PACKAGE="FlashR")
	if (!is.null(vec))
		new("fmV", pointer=vec$pointer, name=vec$name, len=vec$nrow, type=vec$type,
			ele_type=vec$ele_type)
	else
		NULL
}

#' The Normal Distribution
#'
#' \code{fm.rnorm} creates a FlashR vector with random numbers from
#' normal distribution.
#' \code{fm.rnorm.matrix} creates a FlashR matrix with random numbers from
#' normal distribution.
#'
#' If a user provides \code{name} and \code{in.mem} is \code{TRUE}, the created
#' vector/matrix will be kept on disks persistently. That is, even if a user
#' exits from R, the vector/matrix will still be kept on disks.
#'
#' @param n the number of random numbers to be generated.
#' @param nrow the number of rows in the generated matrix.
#' @param ncol the number of columns in the generated matrix.
#' @param mean the mean of the distribution.
#' @param sd the standard deviation of the distribution.
#' @param in.mem whether the vector is stored in memory.
#' @param name the name of the matrix. It's stored on disks, it's used as
#'             the file name.
#' @name fm.rnorm
#' @author Da Zheng <dzheng5@@jhu.edu>
NULL

#' @rdname fm.rnorm
fm.rnorm <- function(n, mean=0, sd=1, in.mem=TRUE, name="")
{
	n <- floor(as.numeric(n))
	if (n <= 0) {
		print("we can't generate a vector of 0 elements")
		return(NULL)
	}
	vec <- .Call("R_FM_create_randmat", "norm", n, 1,
				 as.logical(in.mem), as.character(name),
				 list(mu=as.double(mean), sigma=as.double(sd)), PACKAGE="FlashR")
	if (!is.null(vec))
		new("fmV", pointer=vec$pointer, name=vec$name, len=vec$nrow, type=vec$type,
			ele_type=vec$ele_type)
	else
		NULL
}

#' @rdname fm.runif
fm.runif.matrix <- function(nrow, ncol, min=0, max=1, in.mem=TRUE, name="")
{
	nrow <- floor(as.numeric(nrow))
	ncol <- floor(as.numeric(ncol))
	if (nrow <= 0 || ncol <= 0)
		stop("we can't generate a matrix with 0 rows or cols")

	mat <- .Call("R_FM_create_randmat", "uniform", nrow, ncol,
				 as.logical(in.mem), as.character(name),
				 list(min=as.double(min), max=as.double(max)), PACKAGE="FlashR")
	.new.fm(mat)
}

#' @rdname fm.rnorm
fm.rnorm.matrix <- function(nrow, ncol, mean=0, sd=1, in.mem=TRUE, name="")
{
	nrow <- floor(as.numeric(nrow))
	ncol <- floor(as.numeric(ncol))
	if (nrow <= 0 || ncol <= 0)
		stop("we can't generate a matrix with 0 rows or cols")

	mat <- .Call("R_FM_create_randmat", "norm", nrow, ncol,
				 as.logical(in.mem), as.character(name),
				 list(mu=as.double(mean), sigma=as.double(sd)), PACKAGE="FlashR")
	.new.fm(mat)
}

#' Create a sparse projection matrix.
#'
#' \code{fm.rsparse.proj} creates a sparse projection matrix stored in memory.
#'
#' The non-zero values in the sparse projection matrix are either 1 or -1.
#' Their values are uniformly randomly chosen to be 1 or -1.
#'
#' @param nrow the number of rows in the generated matrix.
#' @param ncol the number of columns in the generated matrix.
#' @param density the ratio of non-zero entries to the total number of elements.
#' @param name the name of the matrix. It's stored on disks, it's used as
#'             the file name.
#' @return a FlashR matrix.
fm.rsparse.proj <- function(nrow, ncol, density, name="")
{
	nrow <- floor(as.numeric(nrow))
	ncol <- floor(as.numeric(ncol))
	if (nrow <= 0 || ncol <= 0)
		stop("we can't generate a matrix with 0 rows or cols")

	mat <- .Call("R_FM_rand_sparse_proj", nrow, ncol, as.numeric(density),
				 PACKAGE="FlashR")
	.new.fm(mat)
}

#' Vector
#'
#' \code{fm.as.vector} converts an object to a FlashR vector.
#' \code{as.vector} converts a FlashR vector to an R vector.
#' \code{fm.is.vector} test whether the input argument is a FlashR vector.
#'
#' Right now, \code{fm.as.vector} can only convert a matrix with only one row
#' or one column. Otherwise, the function returns NULL.
#'
#' @param obj an object
#' @param x a FlashR vector
#' @return a FlashR vector
#' @name vector
NULL

#' @rdname vector
fm.as.vector <- function(obj)
{
	stopifnot(!is.null(obj))
	if (class(obj) == "fmV")
		obj
	else if (class(obj) == "fm") {
		vec <- .Call("R_FM_as_vector", obj, PACKAGE="FlashR")
		if (!is.null(vec))
			.new.fmV(vec)
		else
			NULL
	}
	else if (is.vector(obj))
		fm.conv.R2FM(obj)
	else if (is.matrix(obj))
		fm.conv.R2FM(as.vector(obj))
	else
		NULL
}

#' @rdname vector
setMethod("as.vector", signature(x = "fmV"), function(x) fm.conv.FM2R(x))
#' @rdname vector
setMethod("as.vector", signature(x = "fm"), function(x) as.vector(fm.conv.FM2R(x)))

#' @rdname vector
fm.is.vector <- function(x)
{
	stopifnot(!is.null(x))
	class(x) == "fmV"
}

#' Matrices
#'
#' \code{fm.matrix} creates a matrix from the given set of values.
#' \code{as.matrix} attempts to turn a FlashR matrix to an R matrix.
#' \code{fm.as.matrix} attempts to turn its argument into a FlashR matrix.
#' \code{fm.is.matrix} indicates whether a FlashR object is a matrix.
#'
#' @param x an R object
#' @param vec an R data vector.
#' @param nrow the desired number of rows.
#' @param ncol the desired number of columns.
#' @param byrow logical. If \code{FALSE} (the default) the matrix is filled by
#'	            columns, otherwise the matrix is filled by rows.
#' @name matrix
NULL

#' @rdname matrix
setMethod("as.matrix", signature(x = "fm"), function(x) fm.conv.FM2R(x))
#' @rdname matrix
setMethod("as.matrix", signature(x = "fmV"), function(x)
		  fm.conv.FM2R(fm.as.matrix(x)))

#' @rdname matrix
fm.is.matrix <- function(fm)
{
	stopifnot(!is.null(fm))
	class(fm) == "fm"
}

#' @rdname matrix
fm.as.matrix <- function(x)
{
	stopifnot(!is.null(x))
	if (class(x) == "fm")
		x
	else if (class(x) == "fmV") {
		# A FlashR vector is actually stored in a dense matrix.
		# We only need to construct the fm object in R.
		new("fm", pointer=x@pointer, name=x@name, nrow=x@len,
			ncol=1, type=x@type, ele_type=x@ele_type)
	}
	else {
		# Let's convert it to a FM object
		ret <- fm.conv.R2FM(x)
		# Then try to convert it to a FM matrix.
		fm.as.matrix(ret)
	}
}

#' Convert the data layout of a FlashR matrix.
#'
#' A FMR matrix can store elements in a row-major or column-major order.
#' By changing the data layout, we can improve efficiency of some matrix
#' operations.
#'
#' @param fm a FlashR matrix
#' @param byrow a logical value to determine the data layout of a FlashR matrix.
#' @return a FlashR matrix
#' @author Da Zheng <dzheng5@@jhu.edu>
fm.conv.layout <- function(fm, byrow=FALSE)
{
	stopifnot(!is.null(fm))
	if (class(fm) != "fm")
		return(NULL)

	ret <- .Call("R_FM_conv_layout", fm, as.logical(byrow), PACKAGE="FlashR")
	if (!is.null(ret))
		.new.fm(ret)
	else
		NULL
}

#' Get the data layout of a FlashR matrix.
#'
#' @param fm a FlashR matrix
#' @return a string that indicates the data layout.
#' @author Da Zheng <dzheng5@@jhu.edu>
fm.get.layout <- function(fm)
{
	stopifnot(!is.null(fm))
	if (class(fm) != "fm")
		return(NULL)
	.Call("R_FM_get_layout", fm, PACKAGE="FlashR")
}

#' Convert a regular R object to a FlashR object.
#'
#' If the R object is a matrix, \code{byrow} determines how data in the generated
#' FlashR object is organized in memory.
#'
#' @param obj a regular R object
#' @param byrow a logical value to determine the data layout of a FlashR matrix.
#' @return a FlashR object. If the input R object has 0 element,
#'         \code{fm.conv.R2FM} returns \code{NULL}.
#' @name fm.conv.R2FM
#' @author Da Zheng <dzheng5@@jhu.edu>
fm.conv.R2FM <- function(obj, byrow=FALSE)
{
	stopifnot(!is.null(obj))
	# This function only deals with vectors and matrices of primitive types.
	stopifnot(is.atomic(obj))

	if (length(obj) == 0)
		return(NULL)

	if (is.vector(obj)) {
		vec <- .Call("R_FM_conv_RVec2FM", obj, PACKAGE="FlashR")
		.new.fmV(vec)
	}
	else if(is.matrix(obj)) {
		m <- .Call("R_FM_conv_RMat2FM", obj, as.logical(byrow), PACKAGE="FlashR")
		.new.fm(m)
	}
	else
		NULL
}

#' Convert a FlashR object to a regular R object
#'
#' @param obj a FlashR object
#' @return a regular R object.
#' @name fm.conv.FM2R
#' @author Da Zheng <dzheng5@@jhu.edu>
fm.conv.FM2R <- function(obj)
{
	stopifnot(!is.null(obj))
	if (class(obj) == "fm" && !fm.is.sparse(obj)) {
		nrow <- dim(obj)[1]
		ncol <- dim(obj)[2]
		if (.typeof.int(obj) == "integer")
			ret <- matrix(vector(mode="integer", nrow * ncol), nrow, ncol)
		else if (.typeof.int(obj) == "double")
			ret <- matrix(vector(mode="double", nrow * ncol), nrow, ncol)
		else if (.typeof.int(obj) == "logical")
			ret <- matrix(vector(mode="logical", nrow * ncol), nrow, ncol)
		else {
			print("can't find the type of a dense matrix")
			return(NULL)
		}
		res <- .Call("R_FM_copy_FM2R", obj, ret, PACKAGE="FlashR")
		if (res) ret else NULL
	}
	else if (class(obj) == "fmV") {
		len <- length(obj)
		if (.typeof.int(obj) == "integer")
			ret <- vector(mode="integer", len)
		else if (.typeof.int(obj) == "double")
			ret <- vector(mode="double", len)
		else if (.typeof.int(obj) == "logical")
			ret <- vector(mode="logical", len)
		else {
			print("can't find the type of a vector")
			return(NULL)
		}
		res <- .Call("R_FM_copy_FM2R", obj, ret, PACKAGE="FlashR")
		if (res) ret else NULL
	}
	else if (class(obj) == "fm" && fm.is.sparse(obj)) {
		print("doesn't support convert a sparse matrix to R object")
		NULL
	}
	else {
		print("It has to be a FlashR object")
		NULL
	}
}

#' @rdname matrix
fm.matrix <- function(vec, nrow, ncol, byrow=FALSE)
{
	stopifnot(!is.null(vec))
	if (is.atomic(vec) && length(vec) == 1)
		m <- .Call("R_FM_create_rep_matrix", vec, as.numeric(nrow), as.numeric(ncol),
				   as.logical(byrow), PACKAGE="FlashR")
	else {
		vec <- fm.as.vector(vec)
		if (is.null(vec))
			return(vec)
		m <- .Call("R_FM_create_rep_matrix", vec, as.numeric(nrow), as.numeric(ncol),
				   as.logical(byrow), PACKAGE="FlashR")
	}
	.new.fm(m)
}

#' The information of a FlashR object
#'
#' These functions provide the basic information of a FlashR object.
#'
#' \code{fm.is.sym} indicates whether a matrix is symmetric.
#'
#' \code{fm.matrix.layout} indicates how data in a matrix is organized.
#
#' \code{fm.is.sparse} indicates whether a matrix is sparse.
#'
#' \code{fm.is.sink} indicates whether a FlashR object is a sink matrix.
#'
#' \code{fm.in.mem} indicates whether a FlashR object is stored in memory.
#'
#' @param fm The FlashR object
#' @return \code{fm.is.sym} and \code{fm.is.sparse} returns boolean constants.
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
fm.is.sink <- function(fm)
{
	stopifnot(!is.null(fm))
	if (class(fm) == "fm" || class(fmV) == "fmV")
		.Call("R_FM_is_sink", fm, PACKAGE="FlashR")
	else
		return(FALSE)
}

.typeof.int <- function(fm)
{
	stopifnot(!is.null(fm))
	stopifnot(class(fm) == "fm" || class(fm) == "fmV")
	fm@ele_type
}

#' @rdname fm.info
fm.in.mem <- function(fm)
{
	stopifnot(!is.null(fm))
	stopifnot(class(fm) == "fm" || class(fm) == "fmV")
	.Call("R_FM_is_inmem", fm, PACKAGE="FlashR")
}

#' Convert a FlashR vector to a FlashR factor vector.
#'
#' @param fm a FlashR vector.
#' @param num.levels The number of levels in the factor vector.
#' @return a FlashR factor vector.
fm.as.factor <- function(fm, num.levels = -1)
{
	stopifnot(!is.null(fm))
	if (class(fm) == "fmFactorV")
		fm
	else if (class(fm) == "fmV") {
		if (typeof(fm) != "integer")
			fm <- as.integer(fm)
		if (num.levels < 0) {
			r <- range(fm)
			num.levels <- r[2] - r[1] + 1
			fm <- fm - r[1]
		}
		new("fmFactorV", num.levels=as.integer(num.levels), pointer=fm@pointer,
			name=fm@name, len=fm@len, type=fm@type, ele_type=fm@ele_type)
	}
	else
		stop("The input argument isn't a vector")
}

#' Matrix multiplication
#'
#' Multiply a sparse/dense matrix with a dense vector/matrix.
#'
#' @param fm A FlashR matrix
#' @param mat A FlashR dense matrix.
#' @return a FlashR vector if the second argument is a vector;
#' a FlashR matrix if the second argument is a matrix.
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
		.new.fmV(o)
	else
		.new.fm(o)
}

#' Matrix inner product
#'
#' It takes two operators and performs inner product on a dense matrix
#" and a dense vector/matrix.
#'
#' @param fm A FlashR matrix
#' @param mat A FlashR dense matrix.
#' @param Fun1 The reference or the name of one of the predefined basic binary
#' operators.
#' @param Fun2 The reference or the name of one of the predefined basic binary
#' operators.
#' @return a FlashR vector if the second argument is a vector;
#' a FlashR matrix if the second argument is a matrix.
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
	else {
		o <- .Call("R_FM_inner_prod_dense", fm, mat, Fun1, Fun2,
				   PACKAGE="FlashR")
		if (class(mat) == "fmV") .new.fmV(o) else .new.fm(o)
	}
}

#' The basic operators supported by FlashR.
#'
#' The basic operators are mainly used by the FlashR functions that
#' accept operators as arguments. Such a function includes \code{fm.mapply},
#' \code{fm.inner.prod}, etc.
#'
#' \code{fm.get.basic.op} gets the predefined basic binary operator specified by a user.
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
#' \code{fm.get.basic.uop} gets the predefined basic unary operator specified by a user.
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
#' \code{fm.init.basic.op} initializes the following basic operators.
#' \itemize{
#' \item{\code{fm.bo.add}}{the predifined basic binary operator for addition.}
#' \item{\code{fm.bo.sub}}{the predifined basic binary operator for subtraction.}
#' \item{\code{fm.bo.mul}}{the predifined basic binary operator for multiplication.}
#' \item{\code{fm.bo.div}}{the predifined basic binary operator for division.}
#' \item{\code{fm.bo.min}}{the predifined basic binary operator for computing minimum.}
#' \item{\code{fm.bo.max}}{the predifined basic binary operator for computing maximum.}
#' \item{\code{fm.bo.pow}}{the predifined basic binary operator for computing exponential.}
#' \item{\code{fm.bo.eq}, \code{fm.bo.neq}, \code{fm.bo.gt}, \code{fm.bo.ge},
#'       \code{fm.bo.lt} and \code{fm.bo.le}}
#' {the predefined basic logical operators to compare two elements: ==, >, >=, <, <=.}
#' \item{\code{fm.buo.neg}}{the predefined basic unary operator for negate.}
#' \item{\code{fm.buo.sqrt}}{the predefined basic unary operator for square root.}
#' \item{\code{fm.buo.abs}}{the predefined basic unary operator for absolute value.}
#' \item{\code{fm.buo.not}}{the predefined logical NOT operator.}
#' \item{\code{fm.buo.ceil}}{the predefined basic unary operator of computing
#'       a ceiling of a number.}
#' \item{\code{fm.buo.floor}}{the predefined basic unary operator of computing
#'       a floor of a number.}
#' \item{\code{fm.buo.log}, \code{fm.buo.log2} and \code{fm.buo.log10}}{
#'       the predefined basic unary operators of computing log with different
#'       bases.}
#' \item{\code{fm.buo.round}}{the predefined basic unary operator of rounding
#'       a value.}
#' \item{\code{fm.buo.as.int}}{the predefined basic unary operator of casting
#'       a numeric value to an integer.}
#' \item{\code{fm.buo.as.numeric}}{the predefined basic unary operator of
#'       casting an integer to a numeric value.}
#' }
#'
#' @param name the name of the basic operator.
#' @return a reference to the specified basic operator.
#' @name fm.basic.op
#' @author Da Zheng <dzheng5@@jhu.edu>
fm.get.basic.op <- function(name)
{
	stopifnot(!is.null(name))
	op <- .Call("R_FM_get_basic_op", as.character(name), PACKAGE="FlashR")
	if (!is.null(op))
		new("fm.bo", info=op$info, name=op$name)
}

#' @name fm.basic.op
fm.get.basic.uop <- function(name)
{
	stopifnot(!is.null(name))
	op <- .Call("R_FM_get_basic_uop", as.character(name), PACKAGE="FlashR")
	if (!is.null(op))
		new("fm.bo", info=op$info, name=op$name)
}

#' @name fm.basic.op
fm.bo.add <- NULL
#' @name fm.basic.op
fm.bo.sub <- NULL
#' @name fm.basic.op
fm.bo.mul <- NULL
#' @name fm.basic.op
fm.bo.div <- NULL
#' @name fm.basic.op
fm.bo.min <- NULL
#' @name fm.basic.op
fm.bo.max <- NULL
#' @name fm.basic.op
fm.bo.pow <- NULL
#' @name fm.basic.op
fm.bo.eq <- NULL
#' @name fm.basic.op
fm.bo.neq <- NULL
#' @name fm.basic.op
fm.bo.gt <- NULL
#' @name fm.basic.op
fm.bo.ge <- NULL
#' @name fm.basic.op
fm.bo.lt <- NULL
#' @name fm.basic.op
fm.bo.le <- NULL
#' @name fm.basic.op
fm.bo.or <- NULL
#' @name fm.basic.op
fm.bo.and <- NULL

#' @name fm.basic.op
fm.bo.count <- NULL
#' @name fm.basic.op
fm.bo.which.max <- NULL
#' @name fm.basic.op
fm.bo.which.min <- NULL
#' @name fm.basic.op
fm.bo.euclidean <- NULL

#' @name fm.basic.op
fm.buo.neg <- NULL
#' @name fm.basic.op
fm.buo.sqrt <- NULL
#' @name fm.basic.op
fm.buo.abs <- NULL
#' @name fm.basic.op
fm.buo.not <- NULL
#' @name fm.basic.op
fm.buo.ceil <- NULL
#' @name fm.basic.op
fm.buo.floor <- NULL
#' @name fm.basic.op
fm.buo.log <- NULL
#' @name fm.basic.op
fm.buo.log2 <- NULL
#' @name fm.basic.op
fm.buo.log10 <- NULL
#' @name fm.basic.op
fm.buo.round <- NULL
#' @name fm.basic.op
fm.buo.as.int <- NULL
#' @name fm.basic.op
fm.buo.as.numeric <- NULL

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

#' Create an aggregate operator
#'
#' This function creates an aggregate operator for aggregation operations
#' on a FlashR object.
#'
#' An Aggregate operator has two parts. \code{agg} computes partial
#' aggregation results and \code{combine} combines the partial aggregation
#' results to compute the final result. Both \code{agg} and \code{combine}
#' are the type of \code{fm.basic.op}.
#'
#' The main reason of using two operators is for parallelization. Each thread
#' computes aggregation on part of the object. Eventually, we need an operator
#' to combine all partial aggregation results.
#'
#' For many aggregation operations, \code{agg} and \code{combine} are the same.
#' For example, in the case of summation, both \code{agg} and \code{combine}
#' are simply \code{+}. In some cases, these two operators can be different.
#' For example, when counting the occurences of unique values, \code{agg}
#' is \code{count} and \code{combine} is \code{+}.
#'
#' @param agg a \code{fm.basic.op} operator that computes partial aggregation
#'            results.
#' @param combine a \code{fm.basic.op} operator that computes the final result.
#' @param name a string indicating the name of the aggregation operator.
#' @return a \code{fm.agg.op} operator.
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

.get.apply.op <- function(name)
{
	# TODO we are using to `name' to identify an apply operator now.
	new("fm.apply.op", info=as.integer(0), name=name)
}

#' Aggregation on a FlashR object.
#'
#' This function accepts a basic operator and perform aggregation on
#' the FlashR object with the basic operator.
#'
#' \code{fm.agg} aggregates over the entire object.
#'
#' \code{fm.agg.lazy} aggregates over the entire object, but it performs
#' aggregation lazily.
#'
#' \code{fm.agg.mat} aggregates on the rows or columns of a matrix. It performs
#' aggregation on the shorter dimension lazily, but on the longer dimension
#' immediately.
#'
#' \code{fm.agg.mat.lazy} aggregates on the rows or columns of a matrix and
#' performs aggregation lazily regardless the dimension.
#'
#' @param fm a FlashR object
#' @param op the reference or the name of a predefined basic operator or
#'           the reference to an aggregation operator returned by
#'           \code{fm.create.agg.op}.
#' @param margin the subscript which the function will be applied over.
#' @return \code{fm.agg} returns a scalar, \code{fm.agg.mat} returns
#'         a FlashR vector, \code{fm.agg.lazy} and \code{fm.agg.mat.lazy}
#'         return a FlashR sink matrix.
#' @name fm.agg
fm.agg <- function(fm, op)
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

#' @name fm.agg
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
	.new.fmV(ret)
}

#' @name fm.agg
fm.agg.mat <- function(fm, margin, op)
{
	stopifnot(!is.null(fm) && !is.null(op))
	stopifnot(class(fm) == "fm")
	if (class(op) == "character")
		op <- fm.get.basic.op(op)
	if (class(op) == "fm.bo")
		op <- fm.create.agg.op(op, op, op@name)
	stopifnot(class(op) == "fm.agg.op")
	ret <- .Call("R_FM_agg_mat", fm, as.integer(margin), op, PACKAGE="FlashR")
	.new.fmV(ret)
}

#' @name fm.agg
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
	.new.fmV(ret)
}

.env.int <- new.env()
.env.int$fm.test.na <- TRUE

fm.set.test.na <- function(val)
{
	.env.int$fm.test.na <- val
}

# If one of the inputs has an NA in an element, the corresponding location
# in `res' is set to NA.
.set.na <- function(in1, in2, res)
{
	# This is a special function that helps to test and set NA on
	# the computation results, so it shouldn't call any functions that
	# try to test and set NA on the result.
	if (is.null(res))
		return(NULL)
	else if (!.env.int$fm.test.na)
		return(res)
	else if (.typeof.int(res) == "logical")
		# is.na always return TRUE or FALSE, we don't need to test and set NA
		# on the result. If we do, we'll get infinite recursive calls.
		ifelse(fm.mapply2(is.na(in1), is.na(in2), fm.bo.or,
						  FALSE), NA, res)
	else if (.typeof.int(res) == "integer")
		ifelse(fm.mapply2(is.na(in1), is.na(in2), fm.bo.or,
						  FALSE), as.integer(NA), res)
	else if (.typeof.int(res) == "double")
		ifelse(fm.mapply2(.is.na.only(in1), .is.na.only(in2),
						  fm.bo.or, FALSE), as.double(NA), res)
	else
		# In this case, we don't do anything with the result.
		res
}

# If the input has an NA in an element, the corresponding location in `res'
# is set to NA.
.set.na1 <- function(input, res)
{
	# This is a special function that helps to test and set NA on
	# the computation results, so it shouldn't call any functions that
	# try to test and set NA on the result.
	if (is.null(res))
		return(NULL)
	else if (!.env.int$fm.test.na)
		return(res)
	else if (.typeof.int(res) == "logical")
		# is.na always return TRUE or FALSE, we don't need to test and set NA
		# on the result. If we do, we'll get infinite recursive calls.
		ifelse(is.na(input), NA, res)
	else if (.typeof.int(res) == "integer")
		ifelse(is.na(input), as.integer(NA), res)
	else if (.typeof.int(res) == "double") {
		ifelse(.is.na.only(input), as.double(NA), res)
	}
	else
		# In this case, we don't do anything with the result.
		res
}

.mapply2.fm <- function(o1, o2, FUN, set.na=TRUE)
{
	if (class(FUN) == "character")
		FUN <- fm.get.basic.op(FUN)
	stopifnot(class(FUN) == "fm.bo")
	stopifnot(class(o1) == "fm" || class(o2) == "fm")
	stopifnot(dim(o1)[2] == dim(o2)[2] && dim(o1)[1] == dim(o2)[1])
	ret <- .Call("R_FM_mapply2", FUN, o1, o2, PACKAGE="FlashR")
	ret <- .new.fm(ret)
	if (set.na)
		ret <- .set.na(o1, o2, ret)
	ret
}

.mapply2.fmV <- function(o1, o2, FUN, set.na=TRUE)
{
	if (class(FUN) == "character")
		FUN <- fm.get.basic.op(FUN)
	stopifnot(class(FUN) == "fm.bo")
	stopifnot(class(o1) == "fmV" || class(o2) == "fmV")
	if (length(o1) == length(o2))
		ret <- .Call("R_FM_mapply2", FUN, o1, o2, PACKAGE="FlashR")
	else if (length(o1) == 1)
		return(.mapply2.ANY.fmV(as.vector(o1), o2, FUN, set.na))
	else if (length(o2) == 1)
		return(.mapply2.fmV.ANY(o1, as.vector(o2), FUN, set.na))

	ret <- .new.fmV(ret)
	if (set.na)
		ret <- .set.na(o1, o2, ret)
	ret
}

.mapply2.fm.m <- function(o1, o2, FUN, set.na=TRUE)
{
	if (class(FUN) == "character")
		FUN <- fm.get.basic.op(FUN)
	stopifnot(class(FUN) == "fm.bo")
	stopifnot(class(o1) == "fm")
	stopifnot(nrow(o1) == nrow(o2))
	stopifnot(ncol(o1) == ncol(o2))
	o2 <- fm.conv.R2FM(o2)
	if (is.null(o2))
		return(NULL)
	ret <- .Call("R_FM_mapply2", FUN, o1, o2, PACKAGE="FlashR")
	ret <- .new.fm(ret)
	if (set.na)
		ret <- .set.na(o1, o2, ret)
	ret
}

.mapply2.m.fm <- function(o1, o2, FUN, set.na=TRUE)
{
	if (class(FUN) == "character")
		FUN <- fm.get.basic.op(FUN)
	stopifnot(class(FUN) == "fm.bo")
	stopifnot(class(o2) == "fm")
	stopifnot(nrow(o1) == nrow(o2))
	stopifnot(ncol(o1) == ncol(o2))
	o1 <- fm.conv.R2FM(o1)
	if (is.null(o1))
		return(NULL)
	ret <- .Call("R_FM_mapply2", FUN, o1, o2, PACKAGE="FlashR")
	ret <- .new.fm(ret)
	if (set.na)
		ret <- .set.na(o1, o2, ret)
	ret
}

.mapply2.fm.ANY <- function(o1, o2, FUN, set.na=TRUE)
{
	if (class(FUN) == "character")
		FUN <- fm.get.basic.op(FUN)
	stopifnot(class(FUN) == "fm.bo")
	stopifnot(class(o1) == "fm")
	stopifnot(is.vector(o2))
	if (length(o2) > 1) {
		o2 <- fm.conv.R2FM(o2)
		if (is.null(o2))
			return(NULL)
		fm.mapply.col(o1, o2, FUN, set.na)
	}
	else {
		ret <- .Call("R_FM_mapply2_AE", FUN, o1, o2, PACKAGE="FlashR")
		ret <- .new.fm(ret)
		if (set.na)
			ret <- .set.na1(o1, ret)
		ret
	}
}

.mapply2.ANY.fm <- function(o1, o2, FUN, set.na=TRUE)
{
	if (class(FUN) == "character")
		FUN <- fm.get.basic.op(FUN)
	stopifnot(class(o2) == "fm")
	stopifnot(class(FUN) == "fm.bo")
	if (is.vector(o1) && length(o1) > 1) {
		print("don't support this operation yet.")
		NULL
	}
	else if (is.vector(o1)) {
		ret <- .Call("R_FM_mapply2_EA", FUN, o1, o2, PACKAGE="FlashR")
		ret <- .new.fm(ret)
		if (set.na)
			ret <- .set.na1(o2, ret)
		ret
	}
	else {
		stopifnot(is.matrix(o1))
		o1 <- fm.conv.R2FM(o1)
		.mapply2.fm.fm(o1, o2)
	}
}

.mapply2.fmV.ANY <- function(o1, o2, FUN, set.na=TRUE)
{
	if (class(FUN) == "character")
		FUN <- fm.get.basic.op(FUN)
	stopifnot(class(FUN) == "fm.bo")
	stopifnot(class(o1) == "fmV")
	stopifnot(is.vector(o2))
	if (length(o2) > 1) {
		stopifnot(length(o2) == length(o1))
		o2 <- fm.conv.R2FM(o2)
		if (is.null(o2))
			return(NULL)
		ret <- .Call("R_FM_mapply2", FUN, o1, o2, PACKAGE="FlashR")
	}
	else
		ret <- .Call("R_FM_mapply2_AE", FUN, o1, o2, PACKAGE="FlashR")
	ret <- .new.fmV(ret)
	if (set.na)
		ret <- .set.na1(o1, ret)
	ret
}

.mapply2.ANY.fmV <- function(o1, o2, FUN, set.na=TRUE)
{
	if (class(FUN) == "character")
		FUN <- fm.get.basic.op(FUN)
	stopifnot(class(FUN) == "fm.bo")
	stopifnot(class(o2) == "fmV")
	stopifnot(is.vector(o1))
	if (length(o1) == 1)
		ret <- .Call("R_FM_mapply2_EA", FUN, o1, o2, PACKAGE="FlashR")
	else {
		stopifnot(length(o1) == length(o2))
		o1 <- fm.conv.R2FM(o1)
		if (is.null(o1))
			return(NULL)
		ret <- .Call("R_FM_mapply2", FUN, o1, o2, PACKAGE="FlashR")
	}
	ret <- .new.fmV(ret)
	if (set.na)
		ret <- .set.na1(o2, ret)
	ret
}

#' Apply a Function to two FlashR vectors/matrices.
#'
#' \code{fm.mapply2} applies \code{FUN} to the first elements of each
#' vector/matrix, the second elements, the third elements, and so on.
#' Two vectors/matrices should have the same shape. Currently,
#' \code{fm.mapply2} only accepts predefined basic operators returned
#' by \code{fm.get.basic.op}.
#'
#' \code{fm.mapply.row} and \code{fm.mapply.col} applies to a matrix and
#' a vector. \code{fm.mapply.row} applies \code{FUN} element-wise to
#' each row of the matrix in the left argument and the vector in the right
#' argument. \code{fm.mapply.col} applies \code{FUN} element-wise to
#' each column of the matrix in the left argument and the vector in the right
#' argument.
#'
#' @param o1,o2 a FlashR vector/matrix.
#' @param FUN the reference or the name of one of the predefined basic binary
#' operators.
#' @param set.na a logical value indicating whether to set the value in
#'               the output matrix to NA if the value in the corresponding
#'               location of the input matrix is NA.
#' @return a FlashR vector/matrix.
#' @name fm.mapply2
#' @author Da Zheng <dzheng5@@jhu.edu>
NULL

#' @rdname fm.mapply2
fm.mapply.row <- function(o1, o2, FUN, set.na=TRUE)
{
	if (class(o2) != "fmV")
		o2 <- fm.conv.R2FM(o2)
	stopifnot(class(o1) == "fm" && class(o2) == "fmV")
	if (class(FUN) == "character")
		FUN <- fm.get.basic.op(FUN)
	stopifnot(class(FUN) == "fm.bo")
	if (length(o2) == 1)
		return(fm.mapply2(o1, as.vector(o2), FUN, set.na))

	ret <- .Call("R_FM_mapply2_MV", o1, o2, as.integer(1), FUN,
				 PACKAGE="FlashR")
	ret <- .new.fm(ret)
	if (set.na)
		ret <- .set.na1(o1, ret)
	ret
}

#' @rdname fm.mapply2
fm.mapply.col <- function(o1, o2, FUN, set.na=TRUE)
{
	if (class(o2) != "fmV")
		o2 <- fm.conv.R2FM(o2)
	stopifnot(class(o1) == "fm" && class(o2) == "fmV")
	if (class(FUN) == "character")
		FUN <- fm.get.basic.op(FUN)
	stopifnot(class(FUN) == "fm.bo")
	if (length(o2) == 1)
		return(fm.mapply2(o1, as.vector(o2), FUN, set.na))

	ret <- .Call("R_FM_mapply2_MV", o1, o2, as.integer(2), FUN,
				 PACKAGE="FlashR")
	ret <- .new.fm(ret)
	if (set.na)
		ret <- .set.na1(o1, ret)
	ret
}

setGeneric("fm.mapply2", function(o1, o2, FUN, set.na=TRUE) standardGeneric("fm.mapply2"))
#' @rdname fm.mapply2
setMethod("fm.mapply2", signature(o1 = "fm", o2 = "fm"), .mapply2.fm)
#' @rdname fm.mapply2
setMethod("fm.mapply2", signature(o1 = "fmV", o2 = "fmV"), .mapply2.fmV)
#' @rdname fm.mapply2
setMethod("fm.mapply2", signature(o1 = "fm", o2 = "fmV"), fm.mapply.col)
#' @rdname fm.mapply2
setMethod("fm.mapply2", signature(o1 = "fmV", o2 = "fm"),
		  function(o1, o2, FUN, set.na) {
			  print("This isn't supported currently.")
			  NULL
		  })
#' @rdname fm.mapply2
setMethod("fm.mapply2", signature(o1 = "fm", o2 = "matrix"), .mapply2.fm.m)
#' @rdname fm.mapply2
setMethod("fm.mapply2", signature(o1 = "matrix", o2 = "fm"), .mapply2.m.fm)
#' @rdname fm.mapply2
setMethod("fm.mapply2", signature(o1 = "fm", o2 = "ANY"), .mapply2.fm.ANY)
#' @rdname fm.mapply2
setMethod("fm.mapply2", signature(o1 = "ANY", o2 = "fm"), .mapply2.ANY.fm)
#' @rdname fm.mapply2
setMethod("fm.mapply2", signature(o1 = "fmV", o2 = "ANY"), .mapply2.fmV.ANY)
#' @rdname fm.mapply2
setMethod("fm.mapply2", signature(o1 = "ANY", o2 = "fmV"), .mapply2.ANY.fmV)

.sapply.fm <- function(o, FUN, set.na=TRUE)
{
	stopifnot(class(o) == "fm")
	if (class(FUN) == "character")
		FUN <- fm.get.basic.uop(FUN)
	stopifnot(class(FUN) == "fm.bo")
	ret <- .Call("R_FM_sapply", FUN, o, PACKAGE="FlashR")
	ret <- .new.fm(ret)
	if (as.logical(set.na))
		ret <- .set.na1(o, ret)
	ret
}

.sapply.fmV <- function(o, FUN, set.na=TRUE)
{
	stopifnot(class(o) == "fmV")
	if (class(FUN) == "character")
		FUN <- fm.get.basic.uop(FUN)
	stopifnot(class(FUN) == "fm.bo")
	ret <- .Call("R_FM_sapply", FUN, o, PACKAGE="FlashR")
	ret <- .new.fmV(ret)
	if (as.logical(set.na))
		ret <- .set.na1(o, ret)
	ret
}

#' Apply a Function to a FlashR vector/matrix.
#'
#' \code{sapply} applies \code{FUN} to every element of a vector/matrix.
#' Currently, \code{sapply} only accepts predefined basic operators
#' returned by \code{fm.get.basic.uop}.
#'
#' @param o a FlashR vector/matrix.
#' @param FUN the reference or the name of a predefined uniary operator.
#' @param set.na a logical value indicating whether to set the value in
#'               the output vector or matrix to NA if the value in
#'               the corresponding location of the input vector or matrix
#'               is NA.
#' @return a FlashR vector/matrix.
#' @name fm.sapply
#' @author Da Zheng <dzheng5@@jhu.edu>
setGeneric("fm.sapply", function(o, FUN, set.na=TRUE)  standardGeneric("fm.sapply"))
#' @rdname fm.sapply
setMethod("fm.sapply", signature(o = "fm"), .sapply.fm)
#' @rdname fm.sapply
setMethod("fm.sapply", signature(o = "fmV"), .sapply.fmV)

#' Apply Functions Over Array Margins
#'
#' Apply a predefined function on rows/columns of a FlashR matrix.
#' The predefined function always output a vector of the same length. Thus,
#' the output of \code{fm.apply} is a matrix.
#'
#' Currently, the predefined functions include \code{"rank"} and
#" \code{"sort"}.
#'
#' @param x a FlashR matrix.
#' @param margin an integer. \code{1} indicates rows and \code{2} indicates columns.
#' @param FUN a string that indicates the name of the predefined function.
#' @return a FlashR matrix.
fm.apply <- function(x, margin, FUN)
{
	stopifnot(class(x) == "fm")
	if (class(FUN) == "character")
		FUN <- .get.apply.op(FUN)
	stopifnot(class(FUN) == "fm.apply.bo")
	ret <- .Call("R_FM_apply", FUN, as.integer(margin), x, PACKAGE="FlashR")
	.new.fm(ret)
}

#' Groupby on a FlashR vector.
#'
#' \code{fm.sgroupby} groups elements in a vector based on corresponding
#' \code{labels} and applies \code{FUN} to the elements in each group.
#' \code{FUN} is an aggregation operator.
#'
#' \code{fm.groupby} groups rows/columns of a matrix based on corresponding
#' \code{labels} and applies \code{FUN} to the rows/columns in each group.
#' \code{FUN} is an aggregation operator.
#'
#' @param obj a FlashR vector or matrix
#' @param margin the subscript which the function will be applied over.
#' E.g., for a matrix, \code{1} indicates rows, \code{2} indicates columns.
#' @param factor a FlashR factor vector that indicates how rows/columns
#'               in a matrix should be grouped.
#' @param FUN an aggregation operator returned by \code{fm.create.agg.op}.
#' @return \code{fm.sgroupby} returns a data frame, where the column \code{val}
#' stores all of the unique values in the original data container, and the column
#' \code{agg} stores the aggregate result of the corresponding value.
#' @name fm.groupby
#' @author Da Zheng <dzheng5@@jhu.edu>
NULL

#' @rdname fm.groupby
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
	list(val=.new.fmV(res$val), agg=.new.fmV(res$agg))
}

#' @rdname fm.groupby
fm.groupby <- function(obj, margin, factor, FUN)
{
	stopifnot(class(obj) == "fm")
	if (class(FUN) == "character")
		FUN <- fm.get.basic.op(FUN)
	if (class(FUN) == "fm.bo")
		FUN <- fm.create.agg.op(FUN, FUN, FUN@name)
	stopifnot(class(FUN) == "fm.agg.op")
	stopifnot(class(factor) == "fmFactorV")
	orig.margin <- margin
	if (margin == 1) {
		margin <- 2
		obj <- t(obj)
	}
	res <- .Call("R_FM_groupby", obj, as.integer(margin), factor, FUN,
				 PACKAGE="FlashR")
	if (is.null(res))
		NULL
	else {
		res <- .new.fm(res)
		if (orig.margin == 1) t(res) else res
	}
}

#' Transpose a FlashR matrix.
#'
#' @param m a FlashR matrix
#' @return a FlashR matrix
#' @author Da Zheng <dzheng5@@jhu.edu>
fm.t <- function(m)
{
	stopifnot(!is.null(m))
	stopifnot(class(m) == "fm")
	ret <- .Call("R_FM_transpose", m, PACKAGE="FlashR")
	.new.fm(ret)
}

#' Get a submatrix from a FlashR matrix
#'
#' \code{fm.get.rows} gets specified rows in a FM matrix.
#'
#' \code{fm.get.cols} gets specified columns in a FM matrix.
#'
#' @param fm A FlashR matrix
#' @param idxs an array of column indices in fm.
#' @return a FlashR vector if getting one row or column;
#' a FlashR matrix if getting more than one row or column.
#' @name fm.get.eles
#' @author Da Zheng <dzheng5@@jhu.edu>
fm.get.cols <- function(fm, idxs)
{
	stopifnot(!is.null(fm) && !is.null(idxs))
	stopifnot(class(fm) == "fm")
	if (is.atomic(idxs))
		ret <- .Call("R_FM_get_submat", fm, as.integer(2), as.numeric(idxs),
					 PACKAGE="FlashR")
	else
		ret <- .Call("R_FM_get_submat", fm, as.integer(2), idxs, PACKAGE="FlashR")
	if (!is.null(ret))
		.new.fm(ret)
	else
		NULL
}

#' @rdname fm.get.eles
fm.get.rows <- function(fm, idxs)
{
	stopifnot(!is.null(fm) && !is.null(idxs))
	stopifnot(class(fm) == "fm")
	if (is.atomic(idxs))
		ret <- .Call("R_FM_get_submat", fm, as.integer(1), as.numeric(idxs),
					 PACKAGE="FlashR")
	else
		ret <- .Call("R_FM_get_submat", fm, as.integer(1), idxs, PACKAGE="FlashR")
	if (!is.null(ret))
		.new.fm(ret)
	else
		NULL
}

#' @rdname fm.get.eles
fm.get.eles.vec <- function(fm, idxs)
{
	stopifnot(!is.null(fm) && !is.null(idxs))
	stopifnot(class(fm) == "fmV")
	ret <- .Call("R_FM_get_vec_eles", fm, as.numeric(idxs), PACKAGE="FlashR")
	if (!is.null(ret))
		new("fmV", pointer=ret$pointer, name=ret$name, len=ret$len, type=ret$type,
			ele_type=ret$ele_type)
	else
		NULL
}

#' Materialize virtual FlashR objects.
#'
#' FlashR lazily evaluates many operations and outputs virtual objects
#' that represent computation results. \code{fm.materialize.list} and
#' \code{fm.materialize} explicitly materialize the virtualized computation
#' and save the computation results to memory or disks. Materialization of
#' these virtualized computations triggers materialization of other virtualized
#' computation. By default, FlashR only saves the computation results
#' specified by the arguments of \code{fm.materialize.list} and
#' \code{fm.materialize}.  \code{fm.set.cached} changes the default behavior and
#' notifies FlashR to save the materialized computation results of a virtual
#' matrix in memory or on disks.
#'
#'
#' @param args a list of virtual FlashR objects.
#' @param ... a list of virtual FlashR objects.
#' @param fm a FlashR object.
#' @param in.mem a logical value, indicating whether to save the computation
#'               results in memory.
#' @return a list of materialized compuation results.
#' @name materialize
NULL

#' @rdname materialize
fm.materialize.list <- function(args)
{
	if (length(args) == 0)
		stop("no arguments")
	else if (length(args) == 1) {
		obj <- args[[1]]
		stopifnot(!is.null(obj))
		stopifnot(class(obj) == "fm" || class(obj) == "fmV")
		ret <- .Call("R_FM_materialize", obj, PACKAGE="FlashR")
		if (fm.is.vector(obj))
			.new.fmV(ret)
		else
			.new.fm(ret)
	}
	else {
		for (obj in args) {
			stopifnot(!is.null(obj))
			stopifnot(class(obj) == "fm" || class(obj) == "fmV")
		}
		rets <- .Call("R_FM_materialize_list", args, PACKAGE="FlashR")
		if (is.null(rets))
			return(NULL)
		stopifnot(length(rets) == length(args))
		for (i in 1:length(args)) {
			if (fm.is.vector(args[[i]]))
				rets[[i]] <- .new.fmV(rets[[i]])
			else
				rets[[i]] <- .new.fm(rets[[i]])
		}
		rets
	}
}

#' @rdname materialize
fm.materialize <- function(...)
{
	args <- list(...)
	fm.materialize.list(args)
}

#' @rdname materialize
fm.set.cached <- function(fm, cached, in.mem=fm.in.mem(fm))
{
	stopifnot(!is.null(fm))
	stopifnot(class(fm) == "fm" || class(fm) == "fmV")
	.Call("R_FM_set_materialize_level", fm, as.logical(cached),
		  as.logical(in.mem), PACKAGE="FlashR")
}

#' Write a FlashR object (vector/matrix) to a file
#'
#' @param fm a FlashR object.
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

#' Read a FlashR object (vector/matrix) from a file.
#'
#' @param file a file in the local filesystem.
#' @return a FlashR object (vector/matrix)
#' @author Da Zheng <dzheng5@@jhu.edu>
fm.read.obj <- function(file)
{
	ret <- .Call("R_FM_read_obj", as.character(file), PACKAGE="FlashR")
	if (class(ret) == "fmV")
		.new.fmV(ret)
	else
		.new.fm(ret)
}

#' Convert the Storage of an Object.
#'
#' This function converts the storage of a FlashR vector/matrix.
#' The storage can be memory or disks.
#'
#' If a user provides \code{name} and \code{in.mem} is \code{TRUE},
#' the vector/matrix will be kept on disks persistently. That is, even if a user
#' exits from R, the vector/matrix will still be kept on disks.
#'
#' @param fm a FlashR object.
#' @param in.mem a logical value indicating whether to store the FlashR
#'               object in memory.
#' @param name a string to indicate the name of the new FlashR object.
#' @return a new FlashR object with data stored in specified storage.
fm.conv.store <- function(fm, in.mem, name="")
{
	stopifnot(!is.null(fm))
	stopifnot(class(fm) == "fm" || class(fm) == "fmV")
	ret <- .Call("R_FM_conv_store", fm, as.logical(in.mem),
				 as.character(name), PACKAGE="FlashR")
	if (class(fm) == "fmV")
		.new.fmV(ret)
	else
		.new.fm(ret)
}

#' Combine FlashR Vectors/Matrices by Rows or Columns
#'
#' Take a list of FlashR vectors/matrices and combine them by Columns or rows
#' respectively.
#'
#' @param ... A list of FlashR matrices.
#' @param objs A list of FlashR matrices.
#'
#' All arguments for \code{rbind} and \code{cbind} have to be either all
#' vectors or all matrices.
#'
#' @return A FlashR matrix.
#' @name fm.bind
#' @author Da Zheng <dzheng5@@jhu.edu>
NULL

#' @name fm.bind
fm.rbind.list <- function(objs)
{
	nobjs <- length(objs)
	if (nobjs < 2) {
		return(objs[[1]])
	}
	for (fm in objs) {
		if (class(fm) != "fm" && class(fm) != "fmV") {
			print("fm.rbind only works on FlashR matrix")
			return(NULL)
		}
	}
	ret <- .Call("R_FM_rbind", objs, PACKAGE="FlashR")
	# If it doesn't work, we fall to the default solution and use R
	# to bind matrices.
	if (is.null(ret)) {
		print("Warning! use rbind in R to bind FM matrices")
		for (i in 1:length(objs))
			objs[[i]] <- fm.conv.FM2R(objs[[i]])
		ret <- do.call(rbind, objs)
		fm.conv.R2FM(ret)
	}
	else
		.new.fm(ret)
}

#' @name fm.bind
fm.cbind.list <- function(objs)
{
	nobjs <- length(objs)
	if (nobjs < 2) {
		return(objs[1])
	}
	for (fm in objs) {
		if (class(fm) != "fm" && class(fm) != "fmV") {
			print("fm.cbind only works on FlashR matrix")
			return(NULL)
		}
	}
	ret <- .Call("R_FM_cbind", objs, PACKAGE="FlashR")
	# If it doesn't work, we fall to the default solution and use R
	# to bind matrices.
	if (is.null(ret)) {
		print("Warning! use cbind in R to bind FM matrices")
		for (i in 1:length(objs))
			objs[[i]] <- fm.conv.FM2R(objs[[i]])
		ret <- do.call(cbind, objs)
		fm.conv.R2FM(ret)
	}
	else
		.new.fm(ret)
}

setGeneric("rbind", signature="...")
setGeneric("cbind", signature="...")

#' @name fm.bind
setMethod("rbind", "fm", function(..., deparse.level = 1) {
		  args <- list(...)
		  fm.rbind.list(args)
})
#' @name fm.bind
setMethod("rbind", "fmV", function(..., deparse.level = 1) {
		  args <- list(...)
		  fm.rbind.list(args)
})
#' @name fm.bind
setMethod("cbind", "fm", function(..., deparse.level = 1) {
		  args <- list(...)
		  fm.cbind.list(args)
})
#' @name fm.bind
setMethod("cbind", "fmV", function(..., deparse.level = 1) {
		  args <- list(...)
		  fm.cbind.list(args)
})

#' Conditional Element Selection
#'
#' \code{ifelse} returns a value with the same shape as \code{test} which is
#' filled with elements selected from either \code{yes} or \code{no} depending
#' on whether the element of \code{test} is \code{TRUE} or \code{FALSE}.
#'
#' The current implementation requires either \code{yes} or \code{no} to be
#' a scalar value.
#'
#' @param test a logical FlashR vector or matrix.
#' @param yes a FlashR vector or matrix or an R scalar.
#' @param no a FlashR vector or matrix or an R scalar.
#' @return A FlashR vector or matrix of the same size and attributes
#' (including dimensions) as \code{test} and data values from the values
#' of \code{yes} or \code{no}.
#' @name ifelse
NULL

#' @rdname ifelse
setMethod("ifelse", signature(test = "fm", yes = "fm", no = "ANY"),
		  function(test, yes, no) {
			  ret <- .Call("R_FM_ifelse_no", test, yes, no, PACKAGE="FlashR")
			  .new.fm(ret)
		  })
#' @rdname ifelse
setMethod("ifelse", signature(test = "fm", yes = "ANY", no = "fm"),
		  function(test, yes, no) {
			  ret <- .Call("R_FM_ifelse_yes", test, yes, no, PACKAGE="FlashR")
			  .new.fm(ret)
		  })
#' @rdname ifelse
setMethod("ifelse", signature(test = "fmV", yes = "fmV", no = "ANY"),
		  function(test, yes, no) {
			  ret <- .Call("R_FM_ifelse_no", test, yes, no, PACKAGE="FlashR")
			  .new.fmV(ret)
		  })
#' @rdname ifelse
setMethod("ifelse", signature(test = "fmV", yes = "ANY", no = "fmV"),
		  function(test, yes, no) {
			  ret <- .Call("R_FM_ifelse_yes", test, yes, no, PACKAGE="FlashR")
			  .new.fmV(ret)
		  })

# This returns NA of the right type.
.get.na <- function(type) {
	if (type == "logical")
		NA
	else if (type == "integer")
		as.integer(NA)
	else if (type == "double")
		as.double(NA)
	else {
		stop("unsupported type for NA")
	}
}

.is.na.only <- function(fm)
{
	stopifnot(class(fm) == "fm" || class(fm) == "fmV")
	if (typeof(fm) == "double") {
		ret <- .Call("R_FM_isna", fm, TRUE, PACKAGE="FlashR")
		if (class(fm) == "fm")
			.new.fm(ret)
		else
			.new.fmV(ret)
	}
	else
		is.na(fm)
}

#' 'Not Available' / Missing Values
#'
#' This function indicates which elements are missing.
#'
#' @param x an FlashR vector or matrix.
#' @return a logical FlashR vector or matrix.
#' @name NA
NULL

#' @rdname NA
setMethod("is.na", signature(x = "fm"), function(x) {
		  if (typeof(x) == "double") {
			  ret <- .Call("R_FM_isna", x, FALSE, PACKAGE="FlashR")
			  .new.fm(ret)
		  }
		  else
			  # The result has to be TRUE or FALSE.
			  fm.mapply2(x, .get.na(typeof(x)), fm.bo.eq, FALSE)
		  })
#' @rdname NA
setMethod("is.na", signature(x = "fmV"), function(x) {
		  if (typeof(x) == "double") {
			  ret <- .Call("R_FM_isna", x, FALSE, PACKAGE="FlashR")
			  .new.fmV(ret)
		  }
		  else
			  # The result has to be TRUE or FALSE.
			  fm.mapply2(x, .get.na(typeof(x)), fm.bo.eq, FALSE)
		  })

#' Finite, Infinite and NaN Numbers
#'
#' \code{is.finite} and \code{is.infinite} return a vector of the same length
#' as \code{x}, indicating which elements are finite (not infinite and not
#' missing) or infinite.
#'
#' @param x a FlashR object
#' @return A logical vector of the same length as \code{x}
#' @name is.finite
NULL

#' @rdname is.finite
setMethod("is.nan", signature(x = "fm"), function(x) {
		  if (typeof(x) == "double") {
			  ret <- .Call("R_FM_isnan", x, PACKAGE="FlashR")
			  .new.fm(ret)
		  }
		  else
			  fm.matrix(FALSE, nrow(x), ncol(x))
		  })
#' @rdname is.finite
setMethod("is.nan", signature(x = "fmV"), function(x) {
		  if (typeof(x) == "double") {
			  ret <- .Call("R_FM_isnan", x, PACKAGE="FlashR")
			  .new.fmV(ret)
		  }
		  else
			  fm.rep.int(FALSE, length(x))
		  })

#' @rdname is.finite
setMethod("is.infinite", signature(x = "fm"), function(x) {
		  if (typeof(x) == "double")
			  fm.mapply2(x, Inf, fm.bo.eq, FALSE)
		  else
			  fm.matrix(FALSE, nrow(x), ncol(x))
		  })
#' @rdname is.finite
setMethod("is.infinite", signature(x = "fmV"), function(x) {
		  if (typeof(x) == "double")
			  fm.mapply2(x, Inf, fm.bo.eq, FALSE)
		  else
			  fm.rep.int(FALSE, length(x))
		  })

.test.na.finite <- function(input, res)
{
	if (is.null(res))
		return(NULL)
	else {
		# The result must be float points.
		ifelse(fm.mapply2(is.na(input), is.nan(input), fm.bo.or,
						  FALSE), FALSE, res)
	}
}

#' @rdname is.finite
setMethod("is.finite", signature(x = "fm"), function(x) {
		  if (typeof(x) == "double")
			  .test.na.finite(x, fm.mapply2(x, Inf, fm.bo.neq, FALSE))
		  else
			  fm.matrix(TRUE, nrow(x), ncol(x))
		  })
#' @rdname is.finite
setMethod("is.finite", signature(x = "fmV"), function(x) {
		  if (typeof(x) == "double")
			  .test.na.finite(x, fm.mapply2(x, Inf, fm.bo.neq, FALSE))
		  else
			  fm.rep.int(TRUE, length(x))
		  })

#' The Number of Levels of a Factor
#'
#' Return the number of levels which its argument has.
#'
#' This is applied to a FlashR factor.
#'
#' @param x a FlashR factor
#' @return The length of \code{levels(x)}, which is zero if \code{x}
#' has no levels.
#' @name nlevels
setMethod("nlevels", signature(x = "fmFactorV"), function(x) x@num.levels)

#' Levels Attributes
#'
#' \code{levels} returns the levels attribute of a variable.
#'
#' Currently, levels of a factor vector is a sequence number from 1 to
#' the number of levels. FlashR doesn't support other factor values.
#'
#' @param x a FlashR factor.
#' @return a FlashR vector that contains the value of the levels.
#' @name levels
setMethod("levels", signature(x = "fmFactorV"), function(x)
		  fm.seq.int(1, x@num.levels, 1))

#' Print the information of a FlashR object
#'
#' Print the information of the internal storage representation of
#' a FlashR object. For a virtualized computation result,
#' \code{fm.print.mat.info} prints the internal computation dependancy.
#'
#' @param fm a FlashR object
fm.print.mat.info <- function(fm)
{
	stopifnot(class(fm) == "fm" || class(fm) == "fmV")
	ret <- .Call("R_FM_print_mat_info", fm, PACKAGE="FlashR")
}

#' Google profiler
#'
#' This uses the Google profiler to profile the execution of the code.
#' https://github.com/gperftools/gperftools
#'
#' @param file a string that indicates the file where the profiling result
#'             is saved.
#' @name profile
NULL

#' @rdname profile
fm.start.profiler <- function(file)
{
	if (is.loaded("R_start_profiler"))
		.Call("R_start_profiler", as.character(file), PACKAGE="FlashR")
}

#' @rdname profile
fm.stop.profiler <- function()
{
	if (is.loaded("R_stop_profiler"))
		.Call("R_stop_profiler", PACKAGE="FlashR")
}

.onLoad <- function(libname, pkgname)
{
	library.dynam("FlashR", pkgname, libname, local=FALSE);
	ret <- .Call("R_FM_init", NULL, PACKAGE="FlashR")
	stopifnot(ret)
	fm.init.basic.op()
}
