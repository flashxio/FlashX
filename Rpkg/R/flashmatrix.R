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
#' @param matrix A matrix name.
#' @return A boolean value: true if the matrix has been loaded to FlashMatrixR.
#' @name fm.exist.matrix
#' @author Da Zheng <dzheng5@@jhu.edu>
fm.exist.matrix <- function(fm)
{
	.Call("R_FM_exist_matrix", fm, PACKAGE="FlashGraphR")
}

fm.get.params <- function(name)
{
	.Call("R_FM_get_params", name, PACKAGE="FlashGraphR")
}

#' Load a matrix to FlashMatrixR.
#'
#' Load a matrix to FlashMatrixR from difference sources.
#' 
#' `fm.get.matrix' gets a FlashMatrix object that references a matrix
#' represented by a FlashGraphR object.
#'
#' @param fg A FlashGraphR object.
#' @return a FlashMatrixR object.
#' @name fm.load.matrix
#' @author Da Zheng <dzheng5@@jhu.edu>
#' @examples
#' fg <- fg.load.graph("graph.adj", "graph.index")
#' fm <- fm.get.matrix(fg)
fm.get.matrix <- function(fg)
{
	stopifnot(fg != NULL)
	stopifnot(class(fg) == "fg")
	stopifnot(fg.exist.graph(fg$name))
	ret <- .Call("R_FM_get_matrix_fg", fg, PACKAGE="FlashGraphR")
	structure(ret, class="fm")
}

#' Graph information
#'
#' Functions for providing the basic information of a matrix.
#'
#' `fm.nrow' gets the number of rows in a matrix.
#'
#' `fm.ncol' gets the number of columns in a matrix.
#'
#' `fm.is.sym' indicates whether a matrix is symmetric.
#'
#' @param fm The FlashMatrixR object
#' @return `fm.nrow' and `fm.ncol' returns integer constants.
#' `fm.is.sym' returns boolean constants.
#' @name fm.matrix.info
#' @author Da Zheng <dzheng5@@jhu.edu>

#' @rdname fm.matrix.info
fm.nrow <- function(fm)
{
	stopifnot(fm != NULL)
	stopifnot(class(fm) == "fm")
	fm$nrow
}

#' @rdname fm.graph.info
fm.ncol <- function(fm)
{
	stopifnot(fm != NULL)
	stopifnot(class(fm) == "fm")
	fm$ncol
}

#' @rdname fm.graph.info
fm.is.sym <- function(fm)
{
	stopifnot(fm != NULL)
	stopifnot(class(fm) == "fm")
	fm$sym
}

#' Sparse matrix multiplication
#'
#' Multiply a sparse matrix with a dense vector or a dense matrix.
#'
#' `fm.multiply' multiplies a sparse matrix with a dense vector.
#'
#' `fm.multiply.matrix' multiplies a sparse matrix with a dense matrix.
#'
#' @param fm The FlashMatrixR object
#' @param vec A numueric vector
#' @param m A numeric dense matrix.
#' @return `fm.multiply' returns a numeric vector and `fm.multiply.matrix'
#' returns a numeric dense matrix.
#' @name fm.multiply
#' @author Da Zheng <dzheng5@@jhu.edu>

#' @rdname fm.multiply
fm.multiply <- function(fm, vec)
{
	stopifnot(fm != NULL)
	stopifnot(class(fm) == "fm")
	stopifnot(fm.ncol(fm) == length(vec))
	.Call("R_FM_multiply_v", fm, vec, PACKAGE="FlashGraphR")
}

#' @rdname fm.multiply
fm.multiply.matrix <- function(fm, m, transpose=FALSE)
{
	col.multiply <- function(x) {
		fm.multiply(fm, x, transpose)
	}
	apply(m, 2, col.multiply)
}

print.fm <- function(fm)
{
	stopifnot(fm != NULL)
	stopifnot(class(fm) == "fm")
	directed = "U"
	cat("FlashMatrixR ", fm$name, ": ", fm.nrow(fm), " rows, ", fm.ncol(fm),
		" columns\n", sep="")
}
