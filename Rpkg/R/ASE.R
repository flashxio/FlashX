#' Spectral embedding
#'
#' This computes spectral embedding of a given graph.
#' It supports four variants:
#'
#' Adjacency matrix (A),
#'
#' Augmented matrix (Aug), A - c * D, where c is a constant provided by a user.
#'
#' Laplacian matrix (L), L = D - A,
#'
#' Normalized Laplacian matrix (nL), nL = I - D^(-1/2) * A * D^(-1/2)
#'
#' For the Laplacian matrix and the normalized Laplacian matrix, this function
#' only accepts a symmetric sparse matrix.
#'
#' For the adjacency matrix (A) and augmented matrix (Aug), it computes
#' the largest eigenvalues and the corresponding eigenvectors.
#' For the Laplacian matrix and the normalized Laplacian matrix, it computes
#' the eigenvectors for the smallest eigenvalues.
#'
#' @param fm The FlashMatrixR object
#' @param nev The number of eigenvalues/vectors required.
#' @param which The type of the embedding.
#' @param c The constant used in the Aug variant of embedding.
#' @param ncv The number of vectors in the vector subspace.
#' @param tol Numeric scalar. Stopping criterion: the relative accuracy
#' of the Ritz value is considered acceptable if its error is less than
#' tol times its estimated value. If this is set to zero then machine
#' precision is used.
#' @return A named list with the following members:
#'
#' values
#'
#'    Numeric vector, the desired eigenvalues.
#'
#' vectors
#'
#'    Numeric matrix, the desired eigenvectors as columns. If complex=TRUE
#'   (the default for non-symmetric problems), then the matrix is complex.
#'
#' options
#'
#'   A named list with the supplied options and some information about
#'   the performed calculation, including an ARPACK exit code.
#' @name fm.ase
#' @author Da Zheng <dzheng5@@jhu.edu>
fm.spectral.embedding <- function(fm, nev, which=c("A, Aug, L, nL"),
								  c=1/nrow(fm), ncv=2*nev, tol=1.0e-12)
{
	stopifnot(!is.null(fm))
	stopifnot(class(fm) == "fm")
	# multiply function for eigen on the adjacency matrix
	# this is the default setting.
	multiply <- function(x, extra) fm %*% x
	multiply.right <- function(m) t(fm) %*% m
	multiply.left.diag <- function(v, m) fm.mapply.col(m, v, fm.bo.mul)
	get.degree <- function(fm) fm %*% fm.rep.int(1, ncol(fm))

	directed = !fm.is.sym(fm)
	if (which == "L" || which == "nL") {
		if (directed) {
			print("Don't support embedding on a directed (normalized) Laplacian matrix")
			return(NULL)
		}
	}

	# Compute the opposite of the spectrum.
	comp.oppo <- FALSE
	if (which == "A") {
		if (directed) multiply <- function(x, extra) fm %*% (t(fm) %*% x)
		else multiply <- function(x, extra) fm %*% x
	}
	else if (which == "Aug") {
		cd <- get.degree(fm) * c
		if (directed) {
			multiply <- function(x, extra) {
				t <- t(fm) %*% x + multiply.left.diag(cd, x)
				fm %*% t + multiply.left.diag(cd, t)
			}
			multiply.right <- function(m) t(fm) %*% m + multiply.left.diag(cd, m)
		}
		else
			multiply <- function(x, extra) fm %*% x + multiply.left.diag(cd, x)
	}
	else if (which == "L") {
		d <- get.degree(fm)
		# We compute the largest eigenvalues and then convert them to
		# the smallest eigenvalues. It's easier to compute the largest eigenvalues.
		multiply <- function(x, extra) fm %*% x
		comp.oppo <- TRUE
	}
	else if (which == "nL") {
		d <- 1/sqrt(get.degree(fm))
		# We compute the largest eigenvalues and then convert them to
		# the smallest eigenvalues. It's easier to compute the largest eigenvalues.
		multiply <- function(x, extra) {
			multiply.left.diag(d, fm %*% multiply.left.diag(d, x))
		}
		comp.oppo <- TRUE
	}
	else {
		print("wrong option")
		stopifnot(FALSE)
	}

	ret <- fm.eigen(multiply, sym=TRUE,
					options=list(n=nrow(fm), which="LM",
								 nev=nev, block_size=1, num_blocks=ncv, tol=tol))
	rescale <- function(x) {
		scal <- sqrt(colSums(x * x))
		x <- fm.mapply.row(x, scal, fm.bo.div)
		x <- fm.materialize(x)
	}
	if (directed) {
		ret$vals <- sqrt(ret$vals)
		ret[["left"]] <- ret$vecs
		ret[["right"]] <- rescale(multiply.right(ret$vecs))
		ret$vecs <- NULL
	}
	else if (comp.oppo) {
		if (which == "nL")
			ret$vals <- 1 - ret$vals
		# We can't compute the eigenvalues of the Laplacian matrix.
		else
			ret$vals[1:length(ret$vals)] <- 0
	}
	ret
}
