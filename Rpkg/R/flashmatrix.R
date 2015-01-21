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

#' Create a FlashMatrixR vector
#'
#' `fm.create.vector' creates and initializes a FlashMatrixR vector.
#' It can initialize the vector with a constant `init.v' when `init.rand'
#' is FALSE. Otherwise, it initializes the vector randomly with specified
#' min and max.
#'
#' @param len The length of the vector
#' @param init.v The constat initial value of the vector
#' @param init.rand The boolean indicator whether to initialize the vector
#' randomly.
#' @param init.min The min value of the random generator.
#' @param init.max The max value of the random generator.
#' @return A FlashMatrixR vector
#' @name fm.create.vector
#' @author Da Zheng <dzheng5@@jhu.edu>
fm.create.vector <- function(len, init.v=0, init.rand=FALSE,
							 rand.min=0, rand.max=1)
{
	if (!init.rand)
		ret <- .Call("R_FM_create_vector", as.numeric(len), init.v,
					 PACKAGE="FlashGraphR")
	else {
		# TODO
	}
	structure(ret, class="fmV")
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
	fm$sparse
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
#' @param vec A FlashMatrixR dense vector.
#' @param m A FlashMatrixR matrix
#' @return `fm.multiplyMV' returns a FlashMatrixR vector and
#' `fm.multiplyMM' returns a FlashMatrixR matrix.
#' @name fm.multiply
#' @author Da Zheng <dzheng5@@jhu.edu>

#' @rdname fm.multiply
fm.multiplyMV <- function(fm, vec)
{
	stopifnot(!is.null(fm) && !is.null(vec))
	stopifnot(class(fm) == "fm" && class(vec) == "fmV")
	stopifnot(fm.ncol(fm) == fm.length(vec))
	ret <- .Call("R_FM_multiply_v", fm, vec, PACKAGE="FlashGraphR")
	structure(ret, class="fmV")
}

#' @rdname fm.multiply
fm.multiplyMM <- function(fm, m)
{
	stopifnot(!is.null(fm) && !is.null(m))
	stopifnot(class(fm) == "fm" && class(m) == "fm")
	stopifnot(fm.ncol(fm) == fm.nrow(m))
	ret <- .Call("R_FM_multiply_m", fm, m, PACKAGE="FlashGraphR")
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

print.fm <- function(fm)
{
	stopifnot(!is.null(fm))
	stopifnot(class(fm) == "fm")
	cat("FlashRMatrix ", fm$name, ": ", fm.nrow(fm), " rows, ", fm.ncol(fm),
		" columns, is sparse: ", fm.is.sparse(fm), "\n", sep="")
}

print.fmV <- function(vec)
{
	stopifnot(!is.null(vec))
	stopifnot(class(vec) == "fmV")
	cat("FlashRVector ", vec$name, ": length: ", fm.length(vec), "\n", sep="")
}
