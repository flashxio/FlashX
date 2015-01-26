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
#' @param fg A FlashGraphR object.
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
	ret <- .Call("R_FM_get_matrix_fg", fg, PACKAGE="FlashGraphR")
	structure(ret, class="fm")
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
	ret <- .Call("R_FM_create_vector", as.numeric(times), x,
				 PACKAGE="FlashGraphR")
	structure(ret, class="fmV")
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
	ret <- .Call("R_FM_create_seq", from, to, by,
				 PACKAGE="FlashGraphR")
	structure(ret, class="fmV")
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
	ret <- .Call("R_FM_create_rand", n, min, max,
				 PACKAGE="FlashGraphR")
	structure(ret, class="fmV")
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
		ret <- .Call("R_FM_conv_RVec2FM", obj, PACKAGE="FlashGraphR")
		structure(ret, class="fmV")
	}
	else {
		ret <- .Call("R_FM_conv_RMat2FM", obj, as.logical(byrow), PACKAGE="FlashGraphR")
		structure(ret, class="fm")
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
	if (class(obj) == "fm")
		.Call("R_FM_conv_FM2R", obj, PACKAGE="FlashGraphR")
	else if (class(obj) == "fmV") {
		ret <- .Call("R_FM_conv_FM2R", obj, PACKAGE="FlashGraphR")
		ret[, 1]
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
#' @param data the given set of values.
#' @param nrow the number of rows in the output matrix.
#' @param ncol the number of columns in the output matrix.
#' @param byrow logical. If FALSE (the default) the matrix is filly
#'				columns, otherwise the matrix is filled by rows.
#' @return a FlashMatrixR matrix
fm.matrix <- function(data, nrow, ncol, byrow=FALSE)
{
	stopifnot(!is.null(data))
	stopifnot(class(data) == "fm" || class(data) == "fmV")
	ret <- .Call("R_FM_conv_matrix", data, as.numeric(nrow), as.numeric(ncol),
				 as.logical(byrow), PACKAGE="FlashGraphR")
	structure(ret, class="fm")
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
#
#' `fm.is.sparse' indicates whether a matrix is sparse.
#'
#' `fm.length' gets the length of a FlashMatrixR vector
#'
#' @param fm The FlashMatrixR object
#' @return `fm.nrow' and `fm.ncol' returns numeric constants.
#' `fm.is.sym' and `fm.is.sparse' returns boolean constants.
#' `fm.length' returns a numeric constant.
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
fm.is.sparse <- function(fm)
{
	stopifnot(!is.null(fm))
	stopifnot(class(fm) == "fm")
	fm$type == "sparse"
}

#' @rdname fm.info
fm.length <- function(fm)
{
	stopifnot(!is.null(fm))
	stopifnot(class(fm) == "fmV")
	fm$len
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
	if (class(mat) == "fmV")
		stopifnot(fm.ncol(fm) == fm.length(mat))
	else {
		stopifnot(!fm.is.sparse(mat))
		stopifnot(fm.ncol(fm) == fm.nrow(mat))
	}

	if (fm.is.sparse(fm))
		ret <- .Call("R_FM_multiply_sparse", fm, mat, PACKAGE="FlashGraphR")
	else
		ret <- .Call("R_FM_multiply_dense", fm, mat, PACKAGE="FlashGraphR")

	if (class(mat) == "fmV")
		structure(ret, class="fmV")
	else
		structure(ret, class="fm")
}

#' Sum of a vector/matrix.
#'
#' `sum' returns the sum of all the values in the input vector/matrix.
#'
#' @param m The input vector/matrix
#' @return an integer if the input is an integer vector/matrix; a numeric
#' value if the input is a numeric vector/matrix.
#' @author Da Zheng <dzheng5@@jhu.edu>
fm.sum <- function(m)
{
	stopifnot(!is.null(m))
	stopifnot(class(m) == "fmV" || class(m) == "fm")
	.Call("R_FM_matrix_sum", m, PACKAGE="FlashGraphR")
}

#' The basic operators supported by FlashMatrixR.
#'
#'
#' The basic operators are mainly used by the FlashMatrixR functions that
#' accept operators as arguments. Such a function includes `fm.mapply',
#' `fm.inner.prod', etc.
#'
#' `fm.get.basic.op' gets the predefined basic operator specified by a user.
#' `fm.init.basic.op' defines the basic operators.
#' `fm.bo.add' is the predifined basic operator for addition.
#' `fm.bo.sub' is the predifined basic operator for subtraction.
#' `fm.bo.mul' is the predifined basic operator for multiplication.
#' `fm.bo.div' is the predifined basic operator for division.
#' `fm.bo.min' is the predifined basic operator for computing minimum.
#' `fm.bo.max' is the predifined basic operator for computing maximum.
#' `fm.bo.pow' is the predifined basic operator for computing exponential.
#'
#' @param name the name of the basic operator.
#' @return a reference to the specified basic operator.
#' @name fm.basic.op
#' @author Da Zheng <dzheng5@@jhu.edu>
fm.get.basic.op <- function(name)
{
	stopifnot(!is.null(name))
	ret <- .Call("R_FM_get_basic_op", name, PACKAGE="FlashGraphR")
	structure(ret, class="fm.bo")
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
	class.name <- "";
	if (class(o1) == "fmV") {
		stopifnot(class(o2) == "fmV")
		stopifnot(fm.length(o1) == fm.length(o2))
		class.name <- "fmV"
	}
	else if (class(o1) == "fm") {
		stopifnot(class(o2) == "fm")
		stopifnot(fm.ncol(o1) == fm.ncol(o2) && fm.nrow(o1) == fm.nrow(o2))
		class.name <- "fm"
	}
	else {
		print("o1 has a wrong type")
		return(NULL)
	}

	ret <- .Call("R_FM_mapply2", FUN, o1, o2, PACKAGE="FlashGraphR")
	structure(ret, class=class.name)
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
	ret <- .Call("R_FM_transpose", m, PACKAGE="FlashGraphR")
	structure(ret, class="fm")
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
