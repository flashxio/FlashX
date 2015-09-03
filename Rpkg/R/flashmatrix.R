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
	.Call("R_FM_get_matrix_fg", fg, PACKAGE="FlashGraphR")
}

#' @rdname fm.get.matrix
fm.load.matrix <- function(mat.file, index.file)
{
	.Call("R_FM_load_matrix", mat.file, index.file, PACKAGE="FlashGraphR")
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
	.Call("R_FM_create_vector", as.numeric(times), x, PACKAGE="FlashGraphR")
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
	.Call("R_FM_create_seq", from, to, by, PACKAGE="FlashGraphR")
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
	.Call("R_FM_create_rand", n, min, max, PACKAGE="FlashGraphR")
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
		.Call("R_FM_conv_RVec2FM", obj, PACKAGE="FlashGraphR")
	}
	else {
		.Call("R_FM_conv_RMat2FM", obj, as.logical(byrow), PACKAGE="FlashGraphR")
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
		nrow <- fm.nrow(obj)
		ncol <- fm.ncol(obj)
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
		len <- fm.length(obj)
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
	.Call("R_FM_conv_matrix", vec, as.numeric(nrow), as.numeric(ncol),
				 as.logical(byrow), PACKAGE="FlashGraphR")
}

#' The information of a FlashMatrixR object
#'
#' Functions for providing the basic information of a matrix.
#'
#' `fm.nrow' gets the number of rows in a matrix.
#'
#' `fm.ncol' gets the number of columns in a matrix.
#'
#' `fm.is.sym' indicates whether a matrix is symmetric.
#'
#' `fm.matrix.layout' indicates how data in a matrix is organized.
#
#' `fm.is.sparse' indicates whether a matrix is sparse.
#'
#' `fm.is.vector' indicates whether a FlashMatrixR object is a vector.
#'
#' `fm.length' gets the length of a FlashMatrixR vector
#'
#' `fm.typeof' gets the type of the element in a FlashMatrixR object.
#'
#' @param fm The FlashMatrixR object
#' @return `fm.nrow' and `fm.ncol' returns numeric constants.
#' `fm.is.sym' and `fm.is.sparse' returns boolean constants.
#' `fm.length' returns a numeric constant. `fm.typeof' returns
#' a string.
#' @name fm.info
#' @author Da Zheng <dzheng5@@jhu.edu>

#' @rdname fm.info
fm.nrow <- function(fm)
{
	stopifnot(!is.null(fm))
	stopifnot(class(fm) == "fm")
	fm$nrow
}

#' @rdname fm.info
fm.ncol <- function(fm)
{
	stopifnot(!is.null(fm))
	stopifnot(class(fm) == "fm")
	fm$ncol
}

#' @rdname fm.info
fm.is.sym <- function(fm)
{
	stopifnot(!is.null(fm))
	stopifnot(class(fm) == "fm")
	fm$sym
}

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
	fm$type == "sparse"
}

fm.is.vector <- function(fm)
{
	stopifnot(!is.null(fm))
	class(fm) == "fmV"
}

#' @rdname fm.info
fm.length <- function(fm)
{
	stopifnot(!is.null(fm))
	stopifnot(class(fm) == "fmV")
	fm$len
}

#' @rdname fm.info
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
	.Call("R_FM_as_vector", fm, PACKAGE="FlashGraphR")
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
		stopifnot(fm.ncol(fm) == fm.length(mat))
	}
	else {
		stopifnot(!fm.is.sparse(mat))
		stopifnot(fm.ncol(fm) == fm.nrow(mat))
	}

	if (fm.is.sparse(fm))
		.Call("R_FM_multiply_sparse", fm, mat, PACKAGE="FlashGraphR")
	else
		.Call("R_FM_multiply_dense", fm, mat, PACKAGE="FlashGraphR")
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
#'
#' @param name the name of the basic operator.
#' @return a reference to the specified basic operator.
#' @name fm.basic.op
#' @author Da Zheng <dzheng5@@jhu.edu>
fm.get.basic.op <- function(name)
{
	stopifnot(!is.null(name))
	.Call("R_FM_get_basic_op", name, PACKAGE="FlashGraphR")
}

#' @name fm.basic.op
fm.get.basic.uop <- function(name)
{
	stopifnot(!is.null(name))
	.Call("R_FM_get_basic_uop", name, PACKAGE="FlashGraphR")
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
}

#' Aggregation on a FlashMatrixR object.
#'
#' This function accepts a basic operator and perform aggregation on
#' the FlashMatrixR object with the basic operator.
#'
#' @param bop a basic operator
#' @param fm a FlashMatrixR object
#' @return a scalar
fm.agg <- function(bop, fm)
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
fm.mapply2 <- function(FUN, o1, o2)
{
	stopifnot(!is.null(FUN))
	stopifnot(!is.null(o1) && !is.null(o2))
	stopifnot(class(FUN) == "fm.bo")
	if (class(o1) == "fmV") {
		stopifnot(class(o2) == "fmV")
		stopifnot(fm.length(o1) == fm.length(o2))
	}
	else if (class(o1) == "fm") {
		stopifnot(class(o2) == "fm")
		stopifnot(fm.ncol(o1) == fm.ncol(o2) && fm.nrow(o1) == fm.nrow(o2))
	}
	else {
		print("o1 has a wrong type")
		return(NULL)
	}
	.Call("R_FM_mapply2", FUN, o1, o2, PACKAGE="FlashGraphR")
}

fm.sapply <- function(FUN, o)
{
	stopifnot(!is.null(FUN) && !is.null(o))
	stopifnot(class(FUN) == "fm.bo")
	if (class(o) != "fmV" && class(o) != "fm") {
		print("o has a wrong type")
		return(NULL)
	}
	.Call("R_FM_sapply", FUN, o, PACKAGE="FlashGraphR")
}

fm.ele.wise.op <- function(FUN, o1, o2)
{
	stopifnot(!is.null(FUN))
	stopifnot(!is.null(o1) && !is.null(o2))
	stopifnot(class(FUN) == "fm.bo")

	# This might be the case that o2 is a scalar R variable.
	if (class(o2) != "fmV" && class(o2) != "fm") {
		if (class(o1) != "fmV" && class(o1) != "fm") {
			print("o1 has a wrong type")
			return(NULL)
		}
		# TODO this isn't a very general solution.
		# Can we implement this with sapply?
		.Call("R_FM_mapply2_AE", FUN, o1, o2, PACKAGE="FlashGraphR")
	}
	else if (class(o1) != "fmV" && class(o1) != "fm") {
		if (class(o2) != "fmV" && class(o2) != "fm") {
			print("o2 has a wrong type")
			return(NULL)
		}
		# TODO this isn't a very general solution.
		# Can we implement this with sapply?
		.Call("R_FM_mapply2_EA", FUN, o1, o2, PACKAGE="FlashGraphR")
	}
	else
		fm.mapply2(FUN, o1, o2)
}

#' Transpose a FlashMatrixR matrix.
#'
#' @param m a FlashMatrixR matrix
#' @return a FlashMatrixR matrix
#' @author Da Zheng <dzheng5@@jhu.edu>
fm.t <- function(m)
{
	stopifnot(!is.null(m))
	stopifnot(class(m) == "fm")
	.Call("R_FM_transpose", m, PACKAGE="FlashGraphR")
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
	.Call("R_FM_get_cols", fm, as.integer(idxs), PACKAGE="FlashGraphR")
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
	.Call("R_FM_read_obj", file, PACKAGE="FlashGraphR")
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
	.Call("R_FM_eigen", as.function(func), extra, as.logical(sym),
		  as.list(options), PACKAGE="FlashGraphR")
}

fm.scale <- function(mat, vec, byrow)
{
	.Call("R_FM_scale", mat, vec, as.logical(byrow), PACKAGE="FlashGraphR")
}

print.fm <- function(fm)
{
	stopifnot(!is.null(fm))
	cat("FlashMatrixR matrix ", fm$name, ": ", fm.nrow(fm), " rows, ", fm.ncol(fm),
		" columns, is sparse: ", fm.is.sparse(fm), "\n", sep="")
}

print.fmV <- function(vec)
{
	stopifnot(!is.null(vec))
	cat("FlashVectorR vector ", vec$name, ": length: ", fm.length(vec), "\n", sep="")
}

print.fm.bo <- function(bo)
{
	stopifnot(!is.null(bo))
	cat("FLashMatrixR basic operator:", bo$name, "\n")
}
