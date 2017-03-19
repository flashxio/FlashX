library(FlashR)

context("Test core functions.")

type.set <- c("double", "integer", "logical")
spec.vals <- list(NULL, "NA", "NaN", "Inf", "-Inf")

test_that("load a dense matrix from a text file", {
		  mat.name <- "test.mat"
		  orig.mat <- round(matrix(runif(1000000), 10000, 100) * 100, digits=1)
		  write.table(orig.mat, file="test_mat.csv", sep=",", row.names=FALSE, col.names=FALSE)
		  mat <- fm.load.dense.matrix("test_mat.csv", TRUE, name=mat.name)
		  expect_equal(dim(orig.mat), dim(mat))
		  expect_equal(fm.conv.FM2R(mat), orig.mat)

		  real.ncol <- ncol(mat)
		  mat <- fm.load.dense.matrix("test_mat.csv", TRUE, ncol=real.ncol - 1,
									  name=mat.name)
		  expect_equal(ncol(mat), real.ncol - 1)
		  expect_equal(fm.conv.FM2R(mat), orig.mat[,1:ncol(mat)])
		  mat <- fm.load.dense.matrix("test_mat.csv", TRUE, ncol=real.ncol + 1,
									  name=mat.name)
		  expect_equal(ncol(mat), real.ncol + 1)
		  expect_equal(fm.conv.FM2R(mat)[,1:real.ncol], orig.mat)
		  expect_equal(fm.conv.FM2R(mat[,real.ncol+1]), rep.int(0, nrow(mat)))
		  file.remove("test_mat.csv")
		  # TODO I also need to test when the input text file has missing values.
})

test_that("load a sparse matrix from a text file", {
		  download.file("http://snap.stanford.edu/data/wiki-Vote.txt.gz", "wiki-Vote.txt.gz")
		  system("gunzip wiki-Vote.txt.gz")
		  mat <- fm.load.sparse.matrix("wiki-Vote.txt", in.mem=TRUE, is.sym=FALSE, delim="\t")
		  one <- fm.rep.int(1, nrow(mat))
		  res <- mat %*% one
		  expect_equal(length(res), 8298)
		  expect_equal(sum(res), 103689)

		  mat <- fm.load.sparse.matrix("wiki-Vote.txt", in.mem=FALSE, is.sym=FALSE, delim="\t")
		  one <- fm.rep.int(1, nrow(mat))
		  res <- mat %*% one
		  expect_equal(length(res), 8298)
		  expect_equal(sum(res), 103689)
		  file.remove("wiki-Vote.txt")

		  download.file("http://snap.stanford.edu/data/facebook_combined.txt.gz", "facebook.txt.gz")
		  system("gunzip facebook.txt.gz")
		  mat <- fm.load.sparse.matrix("facebook.txt", in.mem=TRUE, is.sym=TRUE, delim=" ")
		  one <- fm.rep.int(1, nrow(mat))
		  res <- mat %*% one
		  expect_equal(length(res), 4039)
		  expect_equal(sum(res), 176468)

		  mat <- fm.load.sparse.matrix("facebook.txt", in.mem=FALSE, is.sym=TRUE, delim=" ")
		  one <- fm.rep.int(1, nrow(mat))
		  res <- mat %*% one
		  expect_equal(length(res), 4039)
		  expect_equal(sum(res), 176468)
		  file.remove("facebook.txt")

		  download.file("http://snap.stanford.edu/data/soc-LiveJournal1.txt.gz",
						"soc-LiveJournal1.txt.gz")
		  system("gunzip soc-LiveJournal1.txt.gz")
		  mat <- fm.load.sparse.matrix("soc-LiveJournal1.txt", in.mem=TRUE, is.sym=FALSE, delim="\t")
		  one <- fm.rep.int(1, nrow(mat))
		  res <- mat %*% one
		  expect_equal(length(res), 4847571)
		  expect_equal(sum(res), 68475391)

		  mat <- fm.load.sparse.matrix("soc-LiveJournal1.txt", in.mem=FALSE, is.sym=FALSE, delim="\t")
		  one <- fm.rep.int(1, nrow(mat))
		  res <- mat %*% one
		  expect_equal(length(res), 4847571)
		  expect_equal(sum(res), 68475391)
		  file.remove("soc-LiveJournal1.txt")

		  download.file("http://snap.stanford.edu/data/bigdata/communities/com-lj.ungraph.txt.gz",
						"com-lj.ungraph.txt.gz")
		  system("gunzip com-lj.ungraph.txt.gz")
		  mat <- fm.load.sparse.matrix("com-lj.ungraph.txt", in.mem=TRUE, is.sym=TRUE, delim="\t")
		  one <- fm.rep.int(1, nrow(mat))
		  res <- mat %*% one
		  expect_equal(sum(res != 0), 3997962)
		  expect_equal(sum(res), 34681189 * 2)

		  mat <- fm.load.sparse.matrix("com-lj.ungraph.txt", in.mem=FALSE, is.sym=TRUE, delim="\t")
		  one <- fm.rep.int(1, nrow(mat))
		  res <- mat %*% one
		  expect_equal(sum(res != 0), 3997962)
		  expect_equal(sum(res), 34681189 * 2)
		  file.remove("facebook.txt")
})

for (type in type.set) {
test_that(paste("create a vector/matrix with repeat values of", type), {
		  if (type == "double") {
			  vec <- fm.rep.int(as.double(1), 1000000)
			  expect_equal(length(vec), 1000000)
			  expect_equal(typeof(vec), type)
			  expect_true(all(vec == 1))

			  mat <- fm.matrix(as.double(1), 10000, 100)
			  expect_equal(nrow(mat), 10000)
			  expect_equal(ncol(mat), 100)
			  expect_equal(typeof(mat), type)
			  expect_true(all(mat == 1))
		  }
		  else if (type == "integer") {
			  vec <- fm.rep.int(as.integer(1), 1000000)
			  expect_equal(length(vec), 1000000)
			  expect_equal(typeof(vec), type)
			  expect_true(all(vec == 1))

			  mat <- fm.matrix(as.integer(1), 10000, 100)
			  expect_equal(nrow(mat), 10000)
			  expect_equal(ncol(mat), 100)
			  expect_equal(typeof(mat), type)
			  expect_true(all(mat == 1))
		  }
		  else {
			  vec <- fm.rep.int(TRUE, 1000000)
			  expect_equal(length(vec), 1000000)
			  expect_equal(typeof(vec), type)
			  expect_true(all(vec == 1))

			  mat <- fm.matrix(TRUE, 10000, 100)
			  expect_equal(nrow(mat), 10000)
			  expect_equal(ncol(mat), 100)
			  expect_equal(typeof(mat), type)
			  expect_true(all(mat == 1))
		  }
})
}

for (type in type.set) {
test_that(paste("create a vector/matrix with sequence numbers of", type), {
		  if (type == "double") {
			  vec <- fm.seq.int(as.double(1), as.double(1000000), as.double(1))
			  expect_equal(length(vec), 1000000)
			  expect_equal(typeof(vec), type)
			  expect_equal(sum(vec), length(vec) * (length(vec) + 1) / 2)

			  mat <- fm.seq.matrix(as.double(1), as.double(1000000), 10000, 100)
			  expect_equal(nrow(mat), 10000)
			  expect_equal(ncol(mat), 100)
			  expect_equal(typeof(mat), type)
			  expect_equal(sum(mat), length(mat) * (length(mat) + 1) / 2)
		  }
		  else if (type == "integer") {
			  vec <- fm.seq.int(as.integer(1), as.integer(1000000), as.integer(1))
			  expect_equal(length(vec), 1000000)
			  expect_equal(typeof(vec), type)
			  expect_equal(sum(as.double(vec)), length(vec) * (length(vec) + 1) / 2)

			  mat <- fm.seq.matrix(as.integer(1), as.integer(1000000), 10000, 100)
			  expect_equal(nrow(mat), 10000)
			  expect_equal(ncol(mat), 100)
			  expect_equal(typeof(mat), type)
			  expect_equal(sum(as.double(mat)), length(mat) * (length(mat) + 1) / 2)
		  }
})
}

test_that(paste("create a vector/matrix with uniform random numbers"), {
		  vec <- fm.runif(100000000, min=-1, max=1)
		  expect_equal(length(vec), 100000000)
		  expect_equal(typeof(vec), "double")
		  expect_true(mean(vec) < 0.001)

		  mat <- fm.runif.matrix(1000000, 100, min=-1, max=1)
		  expect_equal(nrow(mat), 1000000)
		  expect_equal(ncol(mat), 100)
		  expect_equal(typeof(mat), "double")
		  expect_true(mean(mat) < 0.001)
})

test_that(paste("create a random vector/matrix under normal distribution"), {
		  vec <- fm.rnorm(100000000, mean=0, sd=1)
		  expect_equal(length(vec), 100000000)
		  expect_equal(typeof(vec), "double")
		  expect_true(mean(vec) < 0.001)
		  expect_true(abs(1 - sd(vec)) < 0.001)

		  mat <- fm.rnorm.matrix(1000000, 100, mean=0, sd=1)
		  expect_equal(nrow(mat), 1000000)
		  expect_equal(ncol(mat), 100)
		  expect_equal(typeof(mat), "double")
		  expect_true(mean(vec) < 0.001)
		  expect_true(abs(1 - sd(vec)) < 0.001)
})

data.map <- list()
name.map <- list()

print.depth <- function(depth, ...)
{
	space <- ""
	if (2 - depth > 0) {
		for (i in 1:(2 - depth))
			space <- c(space, "    ")
		space <- paste(space, collapse = '')
	}
	print(paste(space, ...))
}

# This creates in-memory and external-memory vectors of different types.
get.vecs <- function(length, depth, lazy) {
	key <- paste("vec-", length, depth, lazy)
	print.depth(depth, key)
	if (key %in% names(data.map)) {
		print.depth(depth, "The vectors exist")
		return(list(vecs=data.map[[key]], names=name.map[[key]]))
	}

	vecs <- list()
	names <- list()
	if (depth == 0) {
		# Create in-memory vectors.
		logical.vec <- fm.rep.int(TRUE, length)
		int.vec <- fm.seq.int(as.integer(1), as.integer(length), as.integer(1))
		double.vec <- fm.runif(length)
		vecs <- list(logical.vec, int.vec, double.vec)
		names <- list(paste("IM-L-vec-", length, sep=""),
					  paste("IM-I-vec-", length, sep=""),
					  paste("IM-D-vec-", length, sep=""))

		# Create external-memory vectors.
		logical.vec <- fm.conv.store(fm.rep.int(TRUE, length), in.mem=FALSE)
		int.vec <- fm.conv.store(fm.seq.int(as.integer(1), as.integer(length),
											as.integer(1)), in.mem=FALSE)
		double.vec <- fm.runif(length, in.mem=FALSE)
		vecs <- c(vecs, logical.vec, int.vec, double.vec)
		names <- c(names, paste("EM-L-vec-", length, sep=""),
				   paste("EM-I-vec-", length, sep=""),
				   paste("EM-D-vec-", length, sep=""))
	}
	else {
		tmp <- get.mapply.vecs(length, depth, lazy)
		vecs <- c(vecs, tmp$vecs)
		names <- c(names, tmp$names)

		tmp <- get.sapply.vecs(length, depth, lazy)
		vecs <- c(vecs, tmp$vecs)
		names <- c(names, tmp$names)

		tmp <- get.agg.vecs(length, depth, lazy)
		vecs <- c(vecs, tmp$vecs)
		names <- c(names, tmp$names)

		if (!lazy)
			for (i in 1:length(vecs))
				vecs[[i]] <- fm.materialize(vecs[[i]])
	}
	data.map[[key]] <<- vecs
	name.map[[key]] <<- names
	list(vecs=vecs, names=names)
}

# All vectors have the same length.
get.mapply.vecs <- function(length, depth, lazy)
{
	print.depth(depth, "get mapply vec of", length, ", depth:", depth, ", lazy:", lazy)
	res <- list()
	res.names <- list()

	tmp <- get.vecs(length, depth - 1, lazy)
	vecs <- tmp$vecs
	names <- tmp$names
	print.depth(depth, length(vecs), "vecs from the lower level")
	for (i in 1:length(vecs)) {
		for (j in 1:length(vecs)) {
			res <- c(res, vecs[[i]] + vecs[[j]])
			res.names <- c(res.names, paste("+(", names[[i]], ", ",
											names[[j]], ")", sep=""))
		}
	}
	list(vecs=res, names=res.names)
}

get.sapply.vecs <- function(length, depth, lazy)
{
	print.depth(depth, "get sapply vec of", length, ", depth:", depth, ", lazy:", lazy)
	res <- list()
	res.names <- list()
	tmp <- get.vecs(length, depth - 1, lazy)
	vecs <- tmp$vecs
	names <- tmp$names
	print.depth(depth, length(vecs), "vecs from the lower level")
	for (i in 1:length(vecs)) {
		res <- c(res, -vecs[[i]])
		res.names <- c(res.names, paste("-(", names[[i]], ")", sep=""))
	}
	list(vecs=res, names=res.names)
}

get.agg.vecs <- function(length, depth, lazy)
{
	print.depth(depth, "get agg vec of", length, ", depth:", depth, ", lazy:", lazy)
	res <- list()
	res.names <- list()

	tmp <- get.mats(length, 100, depth - 1, lazy)
	mats <- tmp$mats
	names <- tmp$names
	print.depth(depth, length(mats), "mats from the lower level")
	for (i in 1:length(mats)) {
		res <- c(res, fm.as.vector(fm.agg.mat(mats[[i]], 1, "+")))
		res.names <- c(res.names, paste("rowSums(", names[[i]], ")", sep=""))
	}

	if (length < 2000) {
		tmp <- get.mats(length, 100000, depth - 1, lazy)
		mats <- tmp$mats
		names <- tmp$names
		for (i in 1:length(mats)) {
			res <- c(res, fm.as.vector(fm.agg.mat(mats[[i]], 1, "+")))
			res.names <- c(res.names, paste("rowSums(", names[[i]], ")", sep=""))
		}
	}

	tmp <- get.mats(100, length, depth - 1, lazy)
	mats <- tmp$mats
	names <- tmp$names
	for (i in 1:length(mats)) {
		res <- c(res, fm.as.vector(fm.agg.mat(mats[[i]], 2, "+")))
		res.names <- c(res.names, paste("colSums(", names[[i]], ")", sep=""))
	}

	if (length < 2000) {
		tmp <- get.mats(100000, length, depth - 1, lazy)
		mats <- tmp$mats
		names <- tmp$names
		for (i in 1:length(mats)) {
			mat <- fm.agg.mat(mats[[i]], 2, "+")
			res <- c(res, fm.as.vector(mat))
			res.names <- c(res.names, paste("colSums(", names[[i]], ")", sep=""))
		}
	}

	list(vecs=res, names=res.names)
}

IM.large.talls <- list()
IM.large.wides <- list()
EM.large.talls <- list()
EM.large.wides <- list()

# This creates in-memory and external memory of different types.
# It also creates sparse projection matrices.
get.mats <- function(nrow, ncol, depth, lazy)
{
	mats <- list()
	names <- list()

	key <- paste("mat-", nrow, ncol, depth, lazy)
	print.depth(depth, key)
	if (key %in% names(data.map)) {
		print.depth(depth, "The matrices exist")
		return(list(mats=data.map[[key]], names=name.map[[key]]))
	}

	if (depth == 0) {
		# Create in-memory matrices.
		logical.mat <- fm.matrix(TRUE, nrow, ncol)
		int.mat <- fm.seq.matrix(as.integer(1), as.integer(nrow * ncol), nrow, ncol)
		if (nrow > ncol && nrow == 100000) {
			if (ncol %in% names(IM.large.talls))
				double.mat <- IM.large.talls[[ncol + nrow]]
			else {
				double.mat <- fm.runif.matrix(nrow, ncol)
				IM.large.talls[[ncol + nrow]] <- double.mat
			}
		}
		else if (ncol > nrow && ncol == 100000) {
			if (nrow %in% names(IM.large.wides))
				double.mat <- IM.large.wides[[nrow + ncol]]
			else {
				double.mat <- fm.runif.matrix(nrow, ncol)
				IM.large.wides[[nrow + ncol]] <- double.mat
			}
		}
		else
			double.mat <- fm.runif.matrix(nrow, ncol)
		mats <- c(logical.mat, int.mat, double.mat)
		names <- c(paste("IM-L-mat-", nrow, "-", ncol, sep=""),
				   paste("IM-I-mat-", nrow, "-", ncol, sep=""),
				   paste("IM-D-mat-", nrow, "-", ncol, sep=""))

		# Create a sparse project matrix.
		mats <- c(mats, fm.rsparse.proj(nrow, ncol, 0.001))
		names <- c(names, paste("proj-mat-", nrow, "-", ncol, sep=""))

		# Create external-memory matrices
		if (nrow > ncol && nrow == 100000) {
			if (ncol %in% names(EM.large.talls))
				double.mat <- EM.large.talls[[ncol + nrow]]
			else {
				double.mat <- fm.runif.matrix(nrow, ncol, in.mem=FALSE)
				EM.large.talls[[ncol + nrow]] <- double.mat
			}
		}
		else if (ncol > nrow && ncol == 100000) {
			if (nrow %in% names(EM.large.wides))
				double.mat <- EM.large.wides[[nrow + ncol]]
			else {
				double.mat <- fm.runif.matrix(nrow, ncol, in.mem=FALSE)
				EM.large.wides[[nrow + ncol]] <- double.mat
			}
		}
		else
			double.mat <- fm.runif.matrix(nrow, ncol, in.mem=FALSE)
		mats <- c(mats, logical.mat, int.mat, double.mat)
		names <- c(names, paste("EM-L-mat-", nrow, "-", ncol, sep=""),
				   paste("EM-I-mat-", nrow, "-", ncol, sep=""),
				   paste("EM-D-mat-", nrow, "-", ncol, sep=""))
	}
	else {
		tmp <- get.mats(nrow, ncol, depth - 1, lazy)
		mats <- tmp$mats
		names <- tmp$names

		tmp <- get.mapply.mats(nrow, ncol, depth, lazy)
		mats <- c(mats, tmp$mats)
		names <- c(names, tmp$names)

		tmp <- get.sapply.mats(nrow, ncol, depth, lazy)
		mats <- c(mats, tmp$mats)
		names <- c(names, tmp$names)

		tmp <- get.t.mats(nrow, ncol, depth, lazy)
		mats <- c(mats, tmp$mats)
		names <- c(names, tmp$names)

		tmp <- get.groupby.mats(nrow, ncol, depth, lazy)
		mats <- c(mats, tmp$mats)
		names <- c(names, tmp$names)

		tmp <- get.multiply.mats(nrow, ncol, depth, lazy)
		mats <- c(mats, tmp$mats)
		names <- c(names, tmp$names)

		if (!lazy)
			for (i in 1:length(mats))
				mats[[i]] <- fm.materialize(mats[[i]])
		if (!lazy) {
			ret <- list()
			for (mat in mats)
				ret <- c(ret, fm.materialize(mat))
			mats <- ret
		}
	}
	data.map[[key]] <<- mats
	name.map[[key]] <<- names
	list(mats=mats, names=names)
}

# All matrices have the same shape
get.mapply.mats <- function(nrow, ncol, depth, lazy)
{
	print.depth(depth, "get mapply mats of", nrow, "x", ncol,
				", depth:", depth, ", lazy:", lazy)
	tmp <- get.mats(nrow, ncol, depth - 1, lazy)
	mats <- tmp$mats
	names <- tmp$names
	print.depth(depth, "get", length(mats), "mats from the lower level")

	res <- list()
	res.names <- list()
	for (i in 1:length(mats)) {
		for (j in 1:length(mats)) {
			res <- c(res, mats[[i]] + mats[[j]])
			res.names <- c(res.names, paste("+(", names[[i]], ", ",
											names[[j]], ")", sep=""))
		}
	}

	vecs <- get.vecs(nrow(mats[[1]]), depth - 1, lazy)
	for (i in 1:length(mats)) {
		for (j in 1:length(vecs)) {
			res <- c(res, fm.mapply.col(mats[[i]], vecs$vecs[[j]], "+"))
			res.names <- c(res.names, paste("+(", names[[i]], ", ",
											vecs$names[[j]], ")", sep=""))
		}
	}

	vecs <- get.vecs(ncol(mats[[1]]), depth - 1, lazy)
	for (i in 1:length(mats)) {
		for (j in 1:length(vecs)) {
			res <- c(res, fm.mapply.row(mats[[i]], vecs$vecs[[j]], "+"))
			res.names <- c(res.names, paste("+(", names[[i]], ", ",
											vecs$names[[j]], ")", sep=""))
		}
	}

	list(mats=res, names=res.names)
}

get.sapply.mats <- function(nrow, ncol, depth, lazy)
{
	print.depth(depth, "get sapply mats of", nrow, "x", ncol,
				", depth:", depth, ", lazy:", lazy)
	tmp <- get.mats(nrow, ncol, depth - 1, lazy)
	mats <- tmp$mats
	names <- tmp$names
	print.depth(depth, "get", length(mats), "mats from the lower level")

	res <- list()
	res.names <- list()
	for (i in 1:length(mats)) {
		res <- c(res, -mats[[i]])
		res.names <- c(res.names, paste("-(", names[[i]], ")", sep=""))
	}

	list(mats=res, names=res.names)
}

get.t.mats <- function(nrow, ncol, depth, lazy)
{
	print.depth(depth, "get tmats of", nrow, "x", ncol,
				", depth:", depth, ", lazy:", lazy)
	tmp <- get.mats(ncol, nrow, depth - 1, lazy)
	mats <- tmp$mats
	names <- tmp$names
	print.depth(depth, "get", length(mats), "mats from the lower level")

	res <- list()
	res.names <- list()
	for (i in 1:length(mats)) {
		res <- c(res, t(mats[[i]]))
		res.names <- c(res.names, paste("t(", names[[i]], ")", sep=""))
	}

	list(mats=res, names=res.names)
}

get.groupby.mats <- function(nrow, ncol, depth, lazy)
{
	res <- list()
	res.names <- list()

	print.depth(depth, "get groupby mats of", nrow, "x", ncol,
				", depth:", depth, ", lazy:", lazy)
	tmp <- get.mats(nrow * 10, ncol, depth - 1, lazy)
	mats <- tmp$mats
	names <- tmp$names
	print.depth(depth, "get", length(mats), "taller mats from the lower level")
	labels <- fm.as.factor(as.integer(fm.runif(nrow(mats[[1]])) * nrow))
	for (i in 1:length(mats)) {
		mat <- fm.groupby(mats[[i]], 2, labels, "+")
		stopifnot(nrow(mat) == nrow)
		stopifnot(ncol(mat) == ncol)
		res <- c(res, mat)
		res.names <- c(res.names, paste("groupbyRowSums(", names[[i]], ")", sep=""))
	}

	tmp <- get.mats(nrow, ncol * 10, depth - 1, lazy)
	mats <- tmp$mats
	names <- tmp$names
	print.depth(depth, "get", length(mats), "wider mats from the lower level")
	labels <- fm.as.factor(as.integer(fm.runif(ncol(mats[[1]])) * ncol))
	for (i in 1:length(mats)) {
		mat <- fm.groupby(mats[[i]], 1, labels, "+")
		stopifnot(nrow(mat) == nrow)
		stopifnot(ncol(mat) == ncol)
		res <- c(res, mat)
		res.names <- c(res.names, paste("groupbyColSums(", names[[i]], ")", sep=""))
	}

	list(mats=res, names=res.names)
}

get.multiply.mats <- function(nrow, ncol, depth, lazy)
{
	res <- list()
	res.names <- list()

	print.depth(depth, "get multiply mats of", nrow, "x", ncol,
				", depth:", depth, ", lazy:", lazy)
	if (nrow < 2000 && ncol < 2000) {
		tmp <- get.mats(nrow, 100000, depth - 1, lazy)
		mats1 <- tmp$mats
		names1 <- tmp$names
		print.depth(depth, "get", length(mats1), "wide mats from the lower level")
		tmp <- get.mats(100000, ncol, depth - 1, lazy)
		mats2 <- tmp$mats
		names2 <- tmp$names
		print.depth(depth, "get", length(mats2), "tall mats from the lower level")
		for (i in 1:length(mats1)) {
			for (j in 1:length(mats2)) {
				res <- c(res, fm.multiply(mats1[[i]], mats2[[j]]))
				res.names <- c(res.names, paste("crossprod(", names1[[i]], ", ",
												names2[[j]], ")", sep=""))
			}
		}
		for (i in 1:length(mats1)) {
			for (j in 1:length(mats2)) {
				res <- c(res, fm.inner.prod(mats1[[i]], mats2[[j]], "*", "+"))
				res.names <- c(res.names, paste("innerprod(t(", names1[[i]], "), ",
												names2[[j]], ")", sep=""))
			}
		}
	}

	tmp <- get.mats(nrow, 100, depth - 1, lazy)
	mats <- tmp$mats
	names <- tmp$names
	print.depth(depth, "get", length(mats), "tall mats from the lower level")
	# The right matrix will be materialized before multiplication.
	# so we don't need to create virtual matrices.
	tmp <- get.mats(100, ncol, 0, lazy)
	smalls <- tmp$mats
	small.names <- tmp$names
	print.depth(depth, "get", length(smalls), "small mats from the lower level")
	for (i in 1:length(mats)) {
		for (j in 1:length(smalls)) {
			res <- c(res, fm.multiply(mats[[i]], smalls[[j]]))
			res.names <- c(res.names, paste("multiply(", names[[i]], ", ",
											small.names[[j]], ")", sep=""))
		}
	}
	for (i in 1:length(mats)) {
		for (j in 1:length(smalls)) {
			res <- c(res, fm.inner.prod(mats[[i]], smalls[[j]], "*", "+"))
			res.names <- c(res.names, paste("inner prod(", names[[i]], ", ",
											small.names[[j]], ")", sep=""))
		}
	}

	list(mats=res, names=res.names)
}

# This creates matrices of different types stored in memory and on disks.
# These matrices have different sizes.
#	nrow < 16000 && ncol < 16000
#	nrow < 32 || ncol < 32
#	nrow > 100 || ncol > 100
#	nrow > 10000 || ncol > 10000
# It also includes
#	t(mat)
#	mapply(mat)
#	t(mapply(mat))
#	sapply(mat)
#	t(sapply(mat))
#	agg(mat)
#	t(agg(mat))
#	groupby(mat)
#	t(groupby(mat))
#	inner.prod(mat)
#	t(inner.prod(mat))
# TODO
#   rbind/cbind
#   get.rows/get.cols.
#   multiply_tall (rectanglar * rectangular)
#   t(multiply_tall (rectanglar * rectangular))
#   multiply_wide
#   t(multiply_wide)
#   get_rows/cols
#   t(get_rows/cols)

for (depth in 0:2) {
test_that("test lazy evals on matrices", {
		  print(depth)
		  print("mat 1000x100")
		  # Dense matrices.
		  # small matrices.
		  tmp <- get.mats(1000, 100, depth, TRUE)
		  objs.lazy <- tmp$mats
		  cat("#lazy objs:", length(objs.lazy), "\n")
		  obj.names <- tmp$names
		  tmp <- get.mats(1000, 100, depth, FALSE)
		  objs <- tmp$mats
		  cat("#objs:", length(objs), "\n")
		  for (i in 1:length(objs)) {
			  print(i)
			  fm.print.mat.info(objs.lazy[[i]])
			  fm.print.mat.info(objs[[i]])
			  obj.lazy <- objs.lazy[[i]]
			  obj <- objs[[i]]
			  expect_true(sum(abs(obj - obj.lazy)) == 0)
		  }
		  data.map <<- list()
		  name.map <<- list()
		  gc()

		  print("mat 100x100")
		  tmp <- get.mats(100, 100, depth, TRUE)
		  objs.lazy <- tmp$mats
		  obj.names <- tmp$names
		  tmp <- get.mats(100, 100, depth, FALSE)
		  objs <- tmp$mats
		  for (i in 1:length(objs)) {
			  print(i)
			  fm.print.mat.info(objs.lazy[[i]])
			  fm.print.mat.info(objs[[i]])
			  obj.lazy <- objs.lazy[[i]]
			  obj <- objs[[i]]
			  expect_true(sum(abs(obj - obj.lazy)) == 0)
		  }
		  data.map <<- list()
		  name.map <<- list()
		  gc()

		  print("mat 1000x1000")
		  tmp <- get.mats(1000, 1000, depth, TRUE)
		  objs.lazy <- tmp$mats
		  obj.names <- tmp$names
		  tmp <- get.mats(1000, 1000, depth, FALSE)
		  objs <- tmp$mats
		  for (i in 1:length(objs)) {
			  print(i)
			  fm.print.mat.info(objs.lazy[[i]])
			  fm.print.mat.info(objs[[i]])
			  obj.lazy <- objs.lazy[[i]]
			  obj <- objs[[i]]
			  expect_true(sum(abs(obj - obj.lazy)) == 0)
		  }
		  data.map <<- list()
		  name.map <<- list()
		  gc()

		  # regular skinny matrix
		  # regular block matrix
		  # wide block matrix
		  for (ncol in c(10, 100, 500)) {
			  print(paste("tall mat of 100000 x", ncol))
			  tmp <- get.mats(100000, ncol, depth, TRUE)
			  objs.lazy <- tmp$mats
			  obj.names <- tmp$names
			  tmp <- get.mats(100000, ncol, depth, FALSE)
			  objs <- tmp$mats
			  for (i in 1:length(objs)) {
				  print(i)
				  fm.print.mat.info(objs.lazy[[i]])
				  fm.print.mat.info(objs[[i]])
				  obj.lazy <- objs.lazy[[i]]
				  obj <- objs[[i]]
				  expect_true(sum(abs(obj - obj.lazy)) == 0)
			  }
			  data.map <<- list()
			  name.map <<- list()
			  gc()
		  }
})
}


test_that("conversion between R and FlashR", {
		  vec <- runif(1000000)
		  fm.vec <- fm.as.vector(vec)
		  expect_true(fm.is.vector(fm.vec))
		  expect_true(!fm.is.matrix(fm.vec))
		  vec1 <- as.vector(fm.vec)
		  expect_equal(vec, vec1)

		  mat <- matrix(runif(1000000), 100000, 10)
		  fm.mat <- fm.as.matrix(mat)
		  expect_true(!fm.is.vector(fm.mat))
		  expect_true(fm.is.matrix(fm.mat))
		  mat1 <- as.matrix(fm.mat)
		  expect_equal(mat, mat1)
})

test_that("matrix multiply", {
		  mat <- matrix(runif(10000000), 1000000, 10)
		  fm.mat <- fm.as.matrix(mat)
		  prod <- crossprod(mat)
		  fm.prod <- crossprod(fm.mat)
		  expect_true(all(abs((fm.prod - prod) / prod) < 1e-14))

		  mat1 <- matrix(runif(200), 10, 20)
		  fm.mat1 <- fm.as.matrix(mat1)
		  prod <- mat %*% mat1
		  fm.prod <- fm.mat %*% fm.mat1
		  expect_true(all(abs((fm.prod - prod) / prod) < 1e-14))

		  mat <- matrix(runif(100000000), 100000, 1000)
		  fm.mat <- fm.as.matrix(mat)
		  prod <- crossprod(mat)
		  fm.prod <- crossprod(fm.mat)
		  expect_true(all(abs((fm.prod - prod) / prod) < 1e-14))

		  mat1 <- matrix(runif(2000000), 1000, 2000)
		  fm.mat1 <- fm.as.matrix(mat1)
		  prod <- mat %*% mat1
		  fm.prod <- fm.mat %*% fm.mat1
		  expect_true(all(abs((fm.prod - prod) / prod) < 1e-14))
})

test_that("inner product", {
		  mat <- matrix(as.integer(runif(10000000, min=-1, max=1) * 10), 1000000, 10)
		  fm.mat <- fm.as.matrix(mat)
		  prod <- crossprod(mat)
		  fm.prod <- fm.inner.prod(t(fm.mat), fm.mat, "*", "+")
		  expect_equal(fm.conv.FM2R(fm.prod), prod)

		  mat1 <- matrix(as.integer(runif(200, min=-1, max=1) * 10), 10, 20)
		  fm.mat1 <- fm.as.matrix(mat1)
		  prod <- mat %*% mat1
		  fm.prod <- fm.inner.prod(fm.mat, fm.mat1, "*", "+")
		  expect_equal(fm.conv.FM2R(fm.prod), prod)

		  mat <- matrix(as.integer(runif(100000000, min=-1, max=1) * 10), 100000, 1000)
		  fm.mat <- fm.as.matrix(mat)
		  prod <- crossprod(mat)
		  fm.prod <- fm.inner.prod(t(fm.mat), fm.mat, "*", "+")
		  expect_equal(fm.conv.FM2R(fm.prod), prod)

		  mat1 <- matrix(as.integer(runif(2000000, min=-1, max=1) * 10), 1000, 2000)
		  fm.mat1 <- fm.as.matrix(mat1)
		  prod <- mat %*% mat1
		  fm.prod <- fm.inner.prod(fm.mat, fm.mat1, "*", "+")
		  expect_equal(fm.conv.FM2R(fm.prod), prod)
})

test_that("aggregation", {
		  mat <- fm.runif.matrix(10000000, 10)
		  rmat <- as.matrix(mat)
		  sum <- fm.agg(mat, "+")
		  rsum <- sum(mat)
		  expect_equal(fm.conv.FM2R(sum), rsum)

		  mat <- t(mat)
		  sum1 <- fm.agg(mat, "+")
		  expect_equal(fm.conv.FM2R(sum1), rsum)

		  mat <- fm.runif.matrix(100000, 1000)
		  rmat <- as.matrix(mat)
		  sum <- fm.agg(mat, "+")
		  rsum <- sum(mat)
		  expect_equal(fm.conv.FM2R(sum), rsum)

		  mat <- t(mat)
		  sum1 <- fm.agg(mat, "+")
		  expect_equal(fm.conv.FM2R(sum1), rsum)
})

test_that("read/write a dense matrix", {
		  fm.mat <- fm.runif.matrix(100, 20)
		  mat <- fm.conv.FM2R(fm.mat)
		  expect_true(fm.write.obj(fm.mat, "test.mat"))
		  fm.mat1 <- fm.read.obj("test.mat")
		  expect_equal(fm.matrix.layout(fm.mat1), "col")
		  expect_equal(mat, fm.conv.FM2R(fm.mat1))
		  file.remove("test.mat")
})

test_that("read/write a transposed dense matrix", {
		  fm.mat <- fm.runif.matrix(100, 20)
		  expect_true(fm.write.obj(t(fm.mat), "test.mat"))
		  fm.mat1 <- fm.read.obj("test.mat")
		  expect_equal(fm.matrix.layout(fm.mat1), "row")
		  expect_equal(fm.conv.FM2R(t(fm.mat)), fm.conv.FM2R(fm.mat1))
		  file.remove("test.mat")
})

test_that("read/write a dense submatrix", {
		  fm.mat <- fm.runif.matrix(100, 20)
		  expect_true(fm.write.obj(fm.get.cols(fm.mat, 3:12), "test.mat"))
		  fm.mat1 <- fm.read.obj("test.mat")
		  expect_equal(fm.matrix.layout(fm.mat1), "col")
		  expect_equal(fm.conv.FM2R(fm.mat[,3:12]), fm.conv.FM2R(fm.mat1))
		  file.remove("test.mat")
})

test_that("read/write a transposed dense submatrix", {
		  fm.mat <- fm.runif.matrix(100, 20)
		  expect_true(fm.write.obj(t(fm.get.cols(fm.mat, 3:12)), "test.mat"))
		  fm.mat1 <- fm.read.obj("test.mat")
		  expect_equal(fm.matrix.layout(fm.mat1), "row")
		  expect_equal(fm.conv.FM2R(t(fm.mat[,3:12])), fm.conv.FM2R(fm.mat1))
		  file.remove("test.mat")
})

test_that("test groupby rows", {
		  m <- fm.runif.matrix(100, 10)
		  v <- floor(fm.runif(100))
		  labels <- fm.as.factor(as.integer(v), max(v) + 1)
		  agg.sum <- fm.create.agg.op(fm.bo.add, fm.bo.add, "sum")
		  res <- fm.groupby(m, 2, labels, agg.sum)
		  res2 <- fm.agg.mat(m, 2, agg.sum)
		  expect_equal(as.vector(res), as.vector(res2))

		  v <- floor(fm.runif(100, min=0, max=2))
		  labels <- fm.as.factor(as.integer(v), max(v) + 1)
		  agg.sum <- fm.create.agg.op(fm.bo.add, fm.bo.add, "sum")
		  res <- fm.groupby(m, 2, labels, agg.sum)
		  rmat <- fm.conv.FM2R(m)
		  rlabels <- fm.conv.FM2R(v)
		  rres <- matrix(nrow=2, ncol=10)
		  rres[1,] <- colSums(rmat[rlabels == 0,, drop=FALSE])
		  rres[2,] <- colSums(rmat[rlabels == 1,, drop=FALSE])
		  expect_equal(fm.conv.FM2R(res), rres)
})

test_that("test groupby cols", {
		  m <- fm.runif.matrix(100, 10)
		  v <- floor(fm.runif(10))
		  labels <- fm.as.factor(as.integer(v), max(v) + 1)
		  agg.sum <- fm.create.agg.op(fm.bo.add, fm.bo.add, "sum")
		  res <- fm.groupby(m, 1, labels, agg.sum)
		  res2 <- fm.agg.mat(m, 1, agg.sum)
		  expect_equal(as.vector(res), as.vector(res2))

		  v <- floor(fm.runif(10, min=0, max=2))
		  labels <- fm.as.factor(as.integer(v), max(v) + 1)
		  agg.sum <- fm.create.agg.op(fm.bo.add, fm.bo.add, "sum")
		  res <- fm.groupby(m, 1, labels, agg.sum)
		  rmat <- fm.conv.FM2R(m)
		  rlabels <- fm.conv.FM2R(v)
		  rres <- matrix(nrow=100, ncol=2)
		  rres[,1] <- rowSums(rmat[,rlabels == 0, drop=FALSE])
		  rres[,2] <- rowSums(rmat[,rlabels == 1, drop=FALSE])
		  expect_equal(fm.conv.FM2R(res), rres)
})


test_that("test mapply2", {
		  # fm op fm
		  m <- fm.runif.matrix(100, 10)
		  m2 <- fm.runif.matrix(100, 10)
		  res <- fm.mapply2(m, m2, "+", TRUE)
		  rres <- fm.conv.FM2R(m) + fm.conv.FM2R(m2)
		  expect_equal(fm.conv.FM2R(res), rres)

		  # fmV op fmV
		  v <- fm.runif(1000)
		  v2 <- fm.runif(1000)
		  res <- fm.mapply2(v, v2, "+", TRUE)
		  rres <- fm.conv.FM2R(v) + fm.conv.FM2R(v2)
		  expect_equal(fm.conv.FM2R(res), rres)

		  # fm op fmV
		  vcol <- fm.runif(nrow(m))
		  res <- fm.mapply2(m, vcol, "+", TRUE)
		  rres <- fm.conv.FM2R(m) + fm.conv.FM2R(vcol)
		  expect_equal(fm.conv.FM2R(res), rres)

		  # fm op matrix
		  mtmp <- matrix(runif(length(m)), nrow(m), ncol(m))
		  res <- fm.mapply2(m, mtmp, "+", TRUE)
		  rres <- fm.conv.FM2R(m) + mtmp
		  expect_equal(fm.conv.FM2R(res), rres)

		  # matrix op fm
		  res <- fm.mapply2(mtmp, m, "+", TRUE)
		  rres <- mtmp + fm.conv.FM2R(m)
		  expect_equal(fm.conv.FM2R(res), rres)

		  # fm op ANY
		  vcol <- runif(nrow(m))
		  res <- fm.mapply2(m, vcol, "+", TRUE)
		  rres <- fm.conv.FM2R(m) + vcol
		  expect_equal(fm.conv.FM2R(res), rres)

		  res <- fm.mapply2(m, 1, "+", TRUE)
		  rres <- fm.conv.FM2R(m) + 1
		  expect_equal(fm.conv.FM2R(res), rres)

		  # ANY op fm
		  res <- fm.mapply2(1, m, "+", TRUE)
		  rres <- 1 + fm.conv.FM2R(m)
		  expect_equal(fm.conv.FM2R(res), rres)

		  # fmV op ANY
		  vtmp <- runif(length(v))
		  res <- fm.mapply2(v, vtmp, "+", TRUE)
		  rres <- fm.conv.FM2R(v) + vtmp
		  expect_equal(fm.conv.FM2R(res), rres)

		  res <- fm.mapply2(v, 1, "+", TRUE)
		  rres <- fm.conv.FM2R(v) + 1
		  expect_equal(fm.conv.FM2R(res), rres)

		  # ANY op fmV
		  res <- fm.mapply2(vtmp, v, "+", TRUE)
		  rres <- vtmp + fm.conv.FM2R(v)
		  expect_equal(fm.conv.FM2R(res), rres)

		  res <- fm.mapply2(1, v, "+", TRUE)
		  rres <- 1 + fm.conv.FM2R(v)
		  expect_equal(fm.conv.FM2R(res), rres)
})
