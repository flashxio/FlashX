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
							  type="character"))
setClass("fmV", representation(pointer = "externalptr", name = "character",
							   len = "numeric", type="character"))
setClass("fm.bo", representation(info = "integer", name = "character"))

new.fm <- function(fm)
{
	if (!is.null(fm))
		new("fm", pointer=fm$pointer, name=fm$name, nrow=fm$nrow, ncol=fm$ncol,
			type=fm$type)
}

new.fmV <- function(fm)
{
	if (!is.null(fm))
		new("fmV", pointer=fm$pointer, name=fm$name, len=fm$len, type=fm$type)
}

#' Reconfigure FlashMatrixR
#'
#' This reconfigures FlashMatrixR with the settings in the configuration file.
#' The configuration file contains a list of key-value pairs. Each line in
#' the file is a key-value pair in the form of "key_name=value".
#' @param conf.file The configuration file.
#' @name fm.set.conf
#' @author Da Zheng <dzheng5@@jhu.edu>
fm.set.conf <- function(conf.file)
{
	ret <- .Call("R_FM_set_conf", conf.file, PACKAGE="FlashGraphR")
	stopifnot(ret);
}

fm.set.log.level <- function(level)
{
	.Call("R_FM_set_log_level", level, PACKAGE="FlashGraphR")
}

#' Indicate whether a matrix has been loaded to FlashMatrixR
#'
#' This function indicates whether a matrix has been loaded to FlashMatrixR.
#' @param name A matrix name.
#' @return A boolean value: true if the matrix has been loaded to FlashMatrixR.
#' @name fm.exist.matrix
#' @author Da Zheng <dzheng5@@jhu.edu>
fm.exist.matrix <- function(name)
{
	.Call("R_FM_exist_matrix", name, PACKAGE="FlashGraphR")
}

fm.get.params <- function(name)
{
	.Call("R_FM_get_params", name, PACKAGE="FlashGraphR")
}

#' Load a sparse matrix to FlashMatrixR.
#'
#' Load a sparse matrix to FlashMatrixR from difference sources.
#' 
#' `fm.get.matrix' gets a FlashMatrixR matrix that references a matrix
#' represented by a FlashGraphR object.
#'
#' `fm.load.matrix' loads a FlashMatrixR matrix from files.
#' The matrix in the file is in the FlashMatrix format.
#'
#' @param fg A FlashGraphR object.
#' @param mat.file The file that stores the sparse matrix.
#' @param index.file The file that stores the index of the sparse matrix.
#' @return a FlashMatrixR matrix.
#' @name fm.get.matrix
#' @author Da Zheng <dzheng5@@jhu.edu>
#' @examples
#' fg <- fg.load.graph("graph.adj", "graph.index")
#' fm <- fm.get.matrix(fg)
fm.get.matrix <- function(fg)
{
	stopifnot(!is.null(fg))
	stopifnot(class(fg) == "fg")
	stopifnot(fg.exist.graph(fg$name))
	m <- .Call("R_FM_get_matrix_fg", fg, PACKAGE="FlashGraphR")
	new.fm(m)
}

#' @rdname fm.get.matrix
fm.load.matrix <- function(mat, index, t.mat=NULL, t.index=NULL, in.mem=TRUE)
{
	if (is.null(t.mat) || is.null(t.index))
		m <- .Call("R_FM_load_matrix_sym", mat, index, in.mem, PACKAGE="FlashGraphR")
	else
		m <- .Call("R_FM_load_matrix_asym", mat, index, t.mat, t.index, in.mem,
			  PACKAGE="FlashGraphR")
	new.fm(m)
}

#' Create a FlashMatrixR vector with replicated elements.
#'
#' @param x A constant initial value
#' @param times The length of the vector to be generated.
#' @return A FlashMatrixR vector
#' @name fm.rep.int
#' @author Da Zheng <dzheng5@@jhu.edu>
fm.rep.int <- function(x, times)
{
	stopifnot(is.vector(x) && is.atomic(x))
	vec <- .Call("R_FM_create_vector", as.numeric(times), x, PACKAGE="FlashGraphR")
	new.fmV(vec)
}

#' Create a FlashMatrixR vector with a sequence of numbers.
#'
#' @param from the starting value of the sequence.
#' @param to the end value of the sequence.
#' @param by number: increment of the sequence.
#' @return a sequence of numbers that belongs to [from, to]
#' @name fm.seq.int
#' @author Da Zheng <dzheng5@@jhu.edu>
fm.seq.int <- function(from, to, by)
{
	vec <- .Call("R_FM_create_seq", from, to, by, PACKAGE="FlashGraphR")
	new.fmV(vec)
}

#' Create a FlashMatrixR vector with uniformly random numbers.
#'
#' @param n the number of random numbers to be generated.
#' @param min lower limits of the distribution.
#' @param max upper limits of the distribution.
#' @name fm.runif
#' @author Da Zheng <dzheng5@@jhu.edu>
fm.runif <- function(n, min=0, max=1)
{
	vec <- .Call("R_FM_create_rand", n, min, max, PACKAGE="FlashGraphR")
	new.fmV(vec)
}

#' Convert a regular R object to a FlashMatrixR object.
#'
#' If the R object is a matrix, `byrow' determines how data in the generated
#' FlashMatrixR object is organized in memory.
#'
#' @param obj a regular R object
#' @param byrow a logical value to determine the data layout of a FlashMatrixR
#' matrix.
#' @return a FlashMatrixR object.
#' @name fm.conv.R2FM
#' @author Da Zheng <dzheng5@@jhu.edu>
fm.conv.R2FM <- function(obj, byrow=FALSE)
{
	stopifnot(!is.null(obj))
	# This function only deals with vectors and matrices of primitive types.
	stopifnot(is.atomic(obj))

	if (is.vector(obj)) {
		vec <- .Call("R_FM_conv_RVec2FM", obj, PACKAGE="FlashGraphR")
		new.fmV(vec)
	}
	else {
		m <- .Call("R_FM_conv_RMat2FM", obj, as.logical(byrow), PACKAGE="FlashGraphR")
		new.fm(m)
	}
}

#' Convert a FlashMatrixR object to a regular R object
#'
#' @param obj a FlashMatrixR object
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
		.Call("R_FM_copy_FM2R", obj, ret, PACKAGE="FlashGraphR")
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
		.Call("R_FM_copy_FM2R", obj, ret, PACKAGE="FlashGraphR")
		ret
	}
	else {
		print("It has to be a FlashMatrixR object")
		NULL
	}
}

#' Create a matrix from the given set of values.
#'
#' The given set has to be a FlashGraphR matrix or FlashGraphR vector.
#'
#' @param vec a given vector.
#' @param nrow the number of rows in the output matrix.
#' @param ncol the number of columns in the output matrix.
#' @param byrow logical. If FALSE (the default) the matrix is filly
#'				columns, otherwise the matrix is filled by rows.
#' @return a FlashMatrixR matrix
fm.matrix <- function(vec, nrow, ncol, byrow=FALSE)
{
	stopifnot(!is.null(vec))
	stopifnot(class(vec) == "fmV")
	m <- .Call("R_FM_conv_matrix", vec, as.numeric(nrow), as.numeric(ncol),
				 as.logical(byrow), PACKAGE="FlashGraphR")
	new.fm(m)
}

#' The information of a FlashMatrixR object
#'
#' Functions for providing the basic information of a matrix.
#'
#' `fm.is.sym' indicates whether a matrix is symmetric.
#'
#' `fm.matrix.layout' indicates how data in a matrix is organized.
#
#' `fm.is.sparse' indicates whether a matrix is sparse.
#'
#' `fm.is.vector' indicates whether a FlashMatrixR object is a vector.
#'
#' @param fm The FlashMatrixR object
#' @return `fm.is.sym' and `fm.is.sparse' returns boolean constants.
#' @name fm.info
#' @author Da Zheng <dzheng5@@jhu.edu>

#' @rdname fm.info
#fm.is.sym <- function(fm)
#{
#	stopifnot(!is.null(fm))
#	stopifnot(class(fm) == "fm")
#	fm@sym
#}

#' @rdname fm.info
fm.matrix.layout <- function(fm)
{
	stopifnot(!is.null(fm))
	stopifnot(class(fm) == "fm")
	.Call("R_FM_matrix_layout", fm, PACKAGE="FlashGraphR")
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

fm.typeof <- function(fm)
{
	stopifnot(!is.null(fm))
	stopifnot(class(fm) == "fm" || class(fm) == "fmV")
	.Call("R_FM_typeof", fm, PACKAGE="FlashGraphR")
}

#' Convert a FlashMatrixR matrix to a FlashMatrixR vector.
#'
#' The matrix must have only one row or one column. Otherwise, the function
#' returns an error.
#'
#' @param fm a FlashMatrixR matrix
#' @return a FlashMatrixR vector
fm.as.vector <- function(fm)
{
	stopifnot(!is.null(fm))
	stopifnot(class(fm) == "fm")
	vec <- .Call("R_FM_as_vector", fm, PACKAGE="FlashGraphR")
	new.fmV(vec)
}

#' Matrix multiplication
#'
#' Multiply a sparse/dense matrix with a dense vector/matrix.
#'
#' @param fm A FlashMatrixR matrix
#' @param mat A FlashMatrixR dense matrix.
#' @return a FlashMatrixR vector if the second argument is a vector;
#' a FlashMatrixR matrix if the second argument is a matrix.
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
		o <- .Call("R_FM_multiply_sparse", fm, mat, PACKAGE="FlashGraphR")
	else
		o <- .Call("R_FM_multiply_dense", fm, mat, PACKAGE="FlashGraphR")
	if (class(mat) == "fmV")
		new.fmV(o)
	else
		new.fm(o)
}

#' The basic operators supported by FlashMatrixR.
#'
#'
#' The basic operators are mainly used by the FlashMatrixR functions that
#' accept operators as arguments. Such a function includes `fm.mapply',
#' `fm.inner.prod', etc.
#'
#' `fm.get.basic.op' gets the predefined basic binary operator specified by a user.
#'
#' `fm.get.basic.uop' gets the predefined basic unary operator specified by a user.
#'
#' `fm.init.basic.op' defines the basic operators.
#'
#' `fm.bo.add' is the predifined basic binary operator for addition.
#' `fm.bo.sub' is the predifined basic binary operator for subtraction.
#' `fm.bo.mul' is the predifined basic binary operator for multiplication.
#' `fm.bo.div' is the predifined basic binary operator for division.
#' `fm.bo.min' is the predifined basic binary operator for computing minimum.
#' `fm.bo.max' is the predifined basic binary operator for computing maximum.
#' `fm.bo.pow' is the predifined basic binary operator for computing exponential.
#' `fm.bo.eq', `fm.bo.gt' and `fm.bo.ge' are the predefined basic
#' logical operators to compare two elements: ==, >, >=.
#'
#' `fm.buo.neg' is the predefined basic unary operator for negate.
#' `fm.buo.sqrt' is the predefined basic unary operator for square root.
#' `fm.buo.abs' is the predefined basic unary operator for absolute value.
#' `fm.buo.not' is the predefined logical NOT Operator.
#' `fm.buo.ceil' is the predefined basic unary Operator of computing a ceiling
#' of a number.
#' `fm.buo.floor' is the predefined basic unary Operator of computing a floor
#' of a number.
#'
#' @param name the name of the basic operator.
#' @return a reference to the specified basic operator.
#' @name fm.basic.op
#' @author Da Zheng <dzheng5@@jhu.edu>
fm.get.basic.op <- function(name)
{
	stopifnot(!is.null(name))
	op <- .Call("R_FM_get_basic_op", name, PACKAGE="FlashGraphR")
	if (!is.null(op))
		new("fm.bo", info=op$info, name=op$name)
}

#' @name fm.basic.op
fm.get.basic.uop <- function(name)
{
	stopifnot(!is.null(name))
	op <- .Call("R_FM_get_basic_uop", name, PACKAGE="FlashGraphR")
	if (!is.null(op))
		new("fm.bo", info=op$info, name=op$name)
}

#' @name fm.basic.op
fm.init.basic.op <- function()
{
	fm.bo.add <<- fm.get.basic.op("add")
	fm.bo.sub <<- fm.get.basic.op("sub")
	fm.bo.mul <<- fm.get.basic.op("mul")
	fm.bo.div <<- fm.get.basic.op("div")
	fm.bo.min <<- fm.get.basic.op("min")
	fm.bo.max <<- fm.get.basic.op("max")
	fm.bo.pow <<- fm.get.basic.op("pow")
	fm.bo.eq <<- fm.get.basic.op("eq")
	fm.bo.gt <<- fm.get.basic.op("gt")
	fm.bo.ge <<- fm.get.basic.op("ge")

	fm.buo.neg <<- fm.get.basic.uop("neg")
	fm.buo.sqrt <<- fm.get.basic.uop("sqrt")
	fm.buo.abs <<- fm.get.basic.uop("abs")
	fm.buo.not <<- fm.get.basic.uop("not")
	fm.buo.ceil <<- fm.get.basic.uop("ceil")
	fm.buo.floor <<- fm.get.basic.uop("floor")
}

#' Aggregation on a FlashMatrixR object.
#'
#' This function accepts a basic operator and perform aggregation on
#' the FlashMatrixR object with the basic operator.
#'
#' @param fm a FlashMatrixR object
#' @param bop a basic operator
#' @return a scalar
fm.agg <- function(fm, bop)
{
	stopifnot(!is.null(fm) && !is.null(bop))
	stopifnot(class(fm) == "fmV" || class(fm) == "fm")
	stopifnot(class(bop) == "fm.bo")
	.Call("R_FM_agg", bop, fm, PACKAGE="FlashGraphR")
}

#' Apply a Function to two FlashMatrixR vectors/matrices.
#'
#' `mapply' applies `FUN' to the first elements of each vector, the second
#' elements, the third elements, and so on. Two vectors/matrices should have
#' the same shape. Currently, `mapply' only accepts predefined basic operators
#' returned by `fm.get.basic.op'.
#'
#' @param FUN the user-defined operators.
#' @param o1, o2 a FlashMatrixR vector/matrix.
#' @return a FlashMatrixR vector/matrix.
#' @author Da Zheng <dzheng5@@jhu.edu>
setGeneric("fm.mapply2", function(o1, o2, FUN) 0)
setMethod("fm.mapply2", signature(o1 = "fm", o2 = "fm", FUN = "fm.bo"),
		  function(o1, o2, FUN) {
			  stopifnot(dim(o1)[2] == dim(o2)[2] && dim(o1)[1] == dim(o2)[1])
			  ret <- .Call("R_FM_mapply2", FUN, o1, o2, PACKAGE="FlashGraphR")
			  new.fm(ret)
		  })
setMethod("fm.mapply2", signature(o1 = "fmV", o2 = "fmV", FUN = "fm.bo"),
		  function(o1, o2, FUN) {
			  stopifnot(length(o1) == length(o2))
			  ret <- .Call("R_FM_mapply2", FUN, o1, o2, PACKAGE="FlashGraphR")
			  new.fmV(ret)
		  })
setMethod("fm.mapply2", signature(o1 = "fmV", o2 = "ANY", FUN = "fm.bo"),
		  function(o1, o2, FUN) {
			  ret <- .Call("R_FM_mapply2_AE", FUN, o1, o2, PACKAGE="FlashGraphR")
			  new.fmV(ret)
		  })
setMethod("fm.mapply2", signature(o1 = "fm", o2 = "ANY", FUN = "fm.bo"),
		  function(o1, o2, FUN) {
			  ret <- .Call("R_FM_mapply2_AE", FUN, o1, o2, PACKAGE="FlashGraphR")
			  new.fm(ret)
		  })
setMethod("fm.mapply2", signature(o1 = "ANY", o2 = "fmV", FUN = "fm.bo"),
		  function(o1, o2, FUN) {
			  ret <- .Call("R_FM_mapply2_EA", FUN, o1, o2, PACKAGE="FlashGraphR")
			  new.fmV(ret)
		  })
setMethod("fm.mapply2", signature(o1 = "ANY", o2 = "fm", FUN = "fm.bo"),
		  function(o1, o2, FUN) {
			  ret <- .Call("R_FM_mapply2_EA", FUN, o1, o2, PACKAGE="FlashGraphR")
			  new.fm(ret)
		  })

#' Apply a Function to a FlashMatrixR vector/matrix.
#'
#' `sapply' applies `FUN' to every element of a vector/matrix.
#' Currently, `sapply' only accepts predefined basic operators
#' returned by `fm.get.basic.uop'.
#'
#' @param FUN the user-defined operators.
#' @param o a FlashMatrixR vector/matrix.
#' @return a FlashMatrixR vector/matrix.
#' @author Da Zheng <dzheng5@@jhu.edu>
setGeneric("fm.sapply", function(o, FUN) 0)
setMethod("fm.sapply", signature(o = "fm", FUN = "fm.bo"),
		  function(o, FUN) {
			  ret <- .Call("R_FM_sapply", FUN, o, PACKAGE="FlashGraphR")
			  new.fm(ret)
		  })
setMethod("fm.sapply", signature(o = "fmV", FUN = "fm.bo"),
		  function(o, FUN) {
			  ret <- .Call("R_FM_sapply", FUN, o, PACKAGE="FlashGraphR")
			  new.fmV(ret)
		  })


#' Transpose a FlashMatrixR matrix.
#'
#' @param m a FlashMatrixR matrix
#' @return a FlashMatrixR matrix
#' @author Da Zheng <dzheng5@@jhu.edu>
fm.t <- function(m)
{
	stopifnot(!is.null(m))
	stopifnot(class(m) == "fm")
	ret <- .Call("R_FM_transpose", m, PACKAGE="FlashGraphR")
	new.fm(ret)
}

#' Set the specified column of a FlashMatrixR matrix.
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
#	.Call("R_FM_set_cols", fm, as.integer(idxs), m2, PACKAGE="FlashGraphR")
}

#' Get a submatrix composed of some columns from a FlashMatrixR matrix
#'
#' This only works on a column-wise matrix. The submatrix is also organized
#' column-wise. The submatrix shares the same data as the input matrix, so
#' there is no memory copy.
#'
#' @param fm A FlashMatrixR matrix
#' @param idxs an array of column indices in fm.
#' @return a submatrix
#' @author Da Zheng <dzheng5@@jhu.edu>
fm.get.cols <- function(fm, idxs)
{
	stopifnot(!is.null(fm) && !is.null(idxs))
	stopifnot(class(fm) == "fm")
	ret <- .Call("R_FM_get_cols", fm, as.integer(idxs), PACKAGE="FlashGraphR")
	new.fm(ret)
}

fm.materialize <- function(fm)
{
	stopifnot(!is.null(fm))
	stopifnot(class(fm) == "fm" || class(fm) == "fmV")
	ret <- .Call("R_FM_materialize", fm, PACKAGE="FlashGraphR")
	new.fm(ret)
}

#' Write a FlashMatrixR object (vector/matrix) to a file
#'
#' @param fm a FlashMatrixR object.
#' @param file a file in the local filesystem.
#' @return a logical value. True if the object is written to a file
#' successfully. Otherwise, FALSE.
#' @author Da Zheng <dzheng5@@jhu.edu>
fm.write.obj <- function(fm, file)
{
	stopifnot(!is.null(fm))
	stopifnot(class(fm) == "fm" || class(fm) == "fmV")
	.Call("R_FM_write_obj", fm, file, PACKAGE="FlashGraphR")
}

#' Read a FlashMatrixR object (vector/matrix) from a file.
#'
#' @param file a file in the local filesystem.
#' @return a FlashMatrixR object (vector/matrix)
#' @author Da Zheng <dzheng5@@jhu.edu>
fm.read.obj <- function(file)
{
	ret <- .Call("R_FM_read_obj", file, PACKAGE="FlashGraphR")
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
#' solver String. The name of the eigensolver to solve the eigenproblems.
#'				Currently, it supports three eigensolvers: KrylovSchur,
#'				Davidson and LOBPCG. KrylovSchur is the default eigensolver.
#' tol Numeric scalar. Stopping criterion: the relative accuracy of
#'				the Ritz value is considered acceptable if its error is less
#'				than `tol' times its estimated value.
#' block_size Numeric scalar. The eigensolvers use a block extension of an
#'				eigensolver algorithm. The block size determines the number
#'				of the vectors that operate together.
#' num_blocks Numeric scalar. The number of blocks to compute eigenpairs.
#' which Specify which eigenvalues/vectors to compute, character
#'              constant with exactly two characters.
#' Possible values for symmetric input matrices:
#' "LA' Compute 'nev' largest (algebraic) eigenvalues.
#' "SA" Compute "nev" smallest (algebraic) eigenvalues.
#' "LM" Compute `nev' largest (in magnitude) eigenvalues.
#' "SM" Compute `nev' smallest (in magnitude) eigenvalues.
#'
#' @return A named list with the following members:
#'         values: Numeric vector, the desired eigenvalues.
#'         vectors: Numeric matrix, the desired eigenvectors as columns.
#' @name fm.eigen
#' @author Da Zheng <dzheng5@@jhu.edu>
fm.eigen <- function(func, extra=NULL, sym=FALSE, options=NULL,
					 env = parent.frame())
{
	ret <- .Call("R_FM_eigen", as.function(func), extra, as.logical(sym),
		  as.list(options), PACKAGE="FlashGraphR")
	ret$vecs <- new.fm(ret$vecs)
	ret
}

fm.scale <- function(mat, vec, byrow)
{
	ret <- .Call("R_FM_scale", mat, vec, as.logical(byrow), PACKAGE="FlashGraphR")
	new.fm(ret)
}
