library(FlashR)

context("Test small matrices")

get.vec <- function(type, len = 2000, spec.val = NULL, percent = 0)
{
	if (is.null(spec.val)) {
		if (type == "integer")
			return(fm.conv.R2FM(1:len))
		else if (type == "double")
			return(fm.runif(len))
		else if (type == "logical")
			return(fm.conv.R2FM(1:len) < 1000)
		else
			return(NULL)
	}
	else if (spec.val == "NA") {
		if (type == "integer") {
			rvec <- 1:len
			rvec[runif(len) * 100 < percent] <- as.integer(NA)
		}
		else if (type == "double") {
			rvec <- runif(len)
			rvec[runif(len) * 100 < percent] <- as.double(NA)
		}
		else if (type == "logical") {
			rvec <- runif(len) < 0.5
			rvec[runif(len) * 100 < percent] <- as.logical(NA)
		}
		else
			return(NULL)
		return(fm.conv.R2FM(rvec))
	}
	else if (spec.val == "NaN") {
		if (type == "double") {
			rvec <- runif(len)
			rvec[runif(len) * 100 < percent] <- NaN
		}
		else
			return(NULL)
		return(fm.conv.R2FM(rvec))
	}
	else if (spec.val == "Inf") {
		if (type == "double") {
			rvec <- runif(len)
			rvec[runif(len) * 100 < percent] <- Inf
		}
		else
			return(NULL)
		return(fm.conv.R2FM(rvec))
	}
	else if (spec.val == "-Inf") {
		if (type == "double") {
			rvec <- runif(len)
			rvec[runif(len) * 100 < percent] <- -Inf
		}
		else
			return(NULL)
		return(fm.conv.R2FM(rvec))
	}
	else
		return(NULL)
}

get.mat <- function(type, nrow=20, ncol=100, spec.val = NULL, percent=0)
{
	len <- nrow * ncol
	if (is.null(spec.val)) {
		if (type == "integer")
			return(as.integer(fm.seq.matrix(1, nrow * ncol, nrow, ncol)))
		else if (type == "double")
			return(fm.runif.matrix(nrow, ncol))
		else if (type == "logical")
			return(fm.seq.matrix(1, nrow * ncol, nrow, ncol) < 1000)
		else
			return(NULL)
	}
	else if (spec.val == "NA") {
		if (type == "integer") {
			rvec <- 1:len
			rvec[runif(len) * 100 < percent] <- as.integer(NA)
		}
		else if (type == "double") {
			rvec <- runif(len)
			rvec[runif(len) * 100 < percent] <- as.double(NA)
		}
		else if (type == "logical") {
			rvec <- runif(len) < 0.5
			rvec[runif(len) * 100 < percent] <- as.logical(NA)
		}
		else
			return(NULL)
		return(fm.conv.R2FM(matrix(rvec, nrow, ncol)))
	}
	else if (spec.val == "NaN") {
		if (type == "double") {
			rvec <- runif(len)
			rvec[runif(len) * 100 < percent] <- NaN
		}
		else
			return(NULL)
		return(fm.conv.R2FM(matrix(rvec, nrow, ncol)))
	}
	else if (spec.val == "Inf") {
		if (type == "double") {
			rvec <- runif(len)
			rvec[runif(len) * 100 < percent] <- Inf
		}
		else
			return(NULL)
		return(fm.conv.R2FM(matrix(rvec, nrow, ncol)))
	}
	else if (spec.val == "-Inf") {
		if (type == "double") {
			rvec <- runif(len)
			rvec[runif(len) * 100 < percent] <- -Inf
		}
		else
			return(NULL)
		return(fm.conv.R2FM(matrix(rvec, nrow, ncol)))
	}
	else
		return(NULL)
}

get.rvec <- function(type, len = 2000)
{
	if (type == "integer")
		return(1:len)
	else if (type == "double")
		return(runif(len))
	else if (type == "logical")
		return((1:len) < 1000)
	else
		return(NULL)
}

get.rmat <- function(type, nrow=20, ncol=100)
{
	if (type == "integer")
		return(matrix(1:(nrow * ncol), nrow, ncol))
	else if (type == "double")
		return(matrix(runif(nrow * ncol), nrow, ncol))
	else if (type == "logical")
		return(matrix(1:(nrow * ncol), nrow, ncol) < 1000)
	else
		stop("wrong left type")
}

type.set <- c("double", "integer", "logical")
spec.vals <- list(NULL, "NA", "NaN", "Inf", "-Inf")

test_that("Test sequence generator",{
		  fm.vec <- fm.seq.int(1, 10000000, 1)
		  vec <- seq.int(1, 10000000, 1)
		  expect_equal(fm.conv.FM2R(fm.vec), vec)
		  expect_equal(typeof(fm.vec), typeof(vec))

		  fm.vec <- fm.seq.int(as.integer(1), as.integer(10000000), as.integer(1))
		  vec <- seq.int(as.integer(1), as.integer(10000000), as.integer(1))
		  expect_equal(fm.conv.FM2R(fm.vec), vec)
		  expect_equal(typeof(fm.vec), typeof(vec))

		  fm.vec <- fm.seq.int(1, as.integer(10000000), as.integer(1))
		  vec <- seq.int(1, as.integer(10000000), as.integer(1))
		  expect_equal(fm.conv.FM2R(fm.vec), vec)
		  expect_equal(typeof(fm.vec), typeof(vec))
})

test_that("Test random generator", {
		  fm.vec <- fm.runif(10000000, -1, 1)
		  vec <- fm.conv.FM2R(fm.vec)
		  expect_equal(sum(vec >= -1), length(fm.vec))
		  expect_equal(sum(vec <= 1), length(fm.vec))
		  expect_equal(typeof(fm.vec), typeof(vec))
})

for (type in type.set) {
test_that("convert a FlashMatrixR vector to R vector", {
		  fm.vec <- get.vec(type)
		  vec <- fm.conv.FM2R(fm.vec)
		  expect_equal(typeof(fm.vec), type)
		  expect_equal(typeof(vec), type)
		  expect_true(is.vector(vec))
		  expect_true(is.atomic(vec))
		  expect_equal(vec, fm.conv.FM2R(fm.vec))
})
}

for (type in type.set) {
test_that("convert a R vector to FlashMatrixR vector", {
		  vec <- get.rvec(type)
		  fm.vec <- fm.conv.R2FM(vec)
		  expect_equal(typeof(vec), type)
		  expect_equal(typeof(fm.vec), type)
		  expect_equal(length(fm.vec), length(vec))
		  expect_equal(vec, fm.conv.FM2R(fm.vec))
})
}

for (type in type.set) {
test_that("create a column-wise FlashMatrixR matrix", {
		  fm.mat <- get.mat(type, nrow=20, ncol=100)
		  expect_equal(typeof(fm.mat), type)
		  expect_equal(dim(fm.mat)[1], 20)
		  expect_equal(dim(fm.mat)[2], 100)
		  mat <- fm.conv.FM2R(fm.mat)
		  expect_equal(nrow(mat), 20)
		  expect_equal(ncol(mat), 100)
})
}

for (type in type.set) {
name <- paste("convert a", type, "R matrix to a FlashMatrixR matrix")
test_that(name, {
		  mat <- get.rmat(type, nrow=20, ncol=100)
		  expect_equal(typeof(mat), type)
		  fm.mat <- fm.conv.R2FM(mat)
		  expect_equal(dim(fm.mat)[1], 20)
		  expect_equal(dim(fm.mat)[2], 100)
		  expect_equal(fm.conv.FM2R(fm.mat), mat)
})
}

test.MV1 <- function(fm.mat, fm.vec, mat, vec)
{
	res <- mat %*% vec
	fm.res <- fm.mat %*% fm.vec
	expect_equal(res, fm.conv.FM2R(fm.res))
	expect_equal(length(fm.res), dim(fm.mat)[1])
	expect_equal(typeof(fm.res), typeof(res))
}

test.MV2 <- function(fm.vec, fm.mat, vec, mat)
{
	res <- vec %*% mat
	fm.res <- fm.vec %*% fm.mat
	expect_equal(res, fm.conv.FM2R(fm.res))
	expect_equal(length(fm.res), dim(fm.mat)[2])
	expect_equal(typeof(fm.res), typeof(res))
}

test.MM1 <- function(fm.left, fm.right, left, right)
{
	res <- left %*% right
	fm.res <- fm.left %*% fm.right
	expect_equal(res, fm.conv.FM2R(fm.res))
	expect_equal(dim(fm.res)[1], dim(fm.left)[1])
	expect_equal(dim(fm.res)[2], dim(fm.right)[2])
	expect_equal(typeof(fm.res), typeof(res))
}

for (left.type in type.set) {
for (right.type in type.set) {
name <- paste("matrix vector multiply with", left.type, right.type)
test_that(name, {
		  mat <- get.rmat(left.type, nrow=20, ncol=100)
		  vec <- get.rvec(right.type, len=100)
		  stopifnot(typeof(mat) == left.type)
		  stopifnot(typeof(vec) == right.type)
		  test.MV1(fm.conv.R2FM(mat), fm.conv.R2FM(vec), mat, vec)
		  test.MV1(fm.conv.R2FM(mat), vec, mat, vec)
		  test.MV2(fm.conv.R2FM(vec), fm.conv.R2FM(t(mat)), vec, t(mat))
		  test.MV2(vec, fm.conv.R2FM(t(mat)), vec, t(mat))
})
}
}

for (left.type in type.set) {
for (right.type in type.set) {
name <- paste("matrix matrix multiply with", left.type, right.type)
test_that(name, {
		  mat1 <- get.rmat(left.type, nrow=20, ncol=100)
		  mat2 <- get.rmat(right.type, nrow=100, ncol=20)
		  stopifnot(typeof(mat1) == left.type)
		  stopifnot(typeof(mat2) == right.type)
		  test.MM1(fm.conv.R2FM(mat1), fm.conv.R2FM(mat2), mat1, mat2)
		  test.MM1(fm.conv.R2FM(mat1), mat2, mat1, mat2)
		  test.MM1(mat1, fm.conv.R2FM(mat2), mat1, mat2)
})
}
}

test_that("column-wise matrix times a vector", {
		  mat <- matrix(1:2000, 20, 100)
		  vec <- 1:ncol(mat)

		  fm.mat <- fm.conv.R2FM(mat)
		  fm.vec <- fm.conv.R2FM(vec)
		  expect_equal(fm.matrix.layout(fm.mat), "col")
		  test.MV1(fm.mat, fm.vec, mat, vec)
})

test.MM.tmp <- function(left.mat, right.mat)
{
	fm.left <- fm.conv.R2FM(left.mat)
	fm.right <- fm.conv.R2FM(right.mat)
	test.MM1(fm.left, fm.right, left.mat, right.mat)
}

test_that("matrix multiply: wide vs. tall", {
		  left.mat <- matrix(runif(2000), 20, 100)
		  right.mat <- matrix(runif(2000), 100, 20)
		  test.MM.tmp(left.mat, right.mat)
})

test_that("matrix multiply: tall vs. small", {
		  left.mat <- matrix(runif(2000), 100, 20)
		  right.mat <- matrix(runif(2000), 20, 10)
		  test.MM.tmp(left.mat, right.mat)
})

# TODO we need to test `^`
bin.ops <- list(`+`, `-`, `*`, `/`, `==`, `!=`, `>`, `>=`, `<`, `<=`, `|`, `&`, `^`,
				pmin, pmax)
bin.op.strs <- list("+", "-", "*", "/", "==", "!=", ">", ">=", "<", "<=", "|", "&", "^",
				"pmin", "pmax")

# Test element-wise vector vector operation.
for (left.spec in spec.vals) {
for (right.spec in spec.vals) {
for (i in 1:length(bin.ops)) {
	bin.op <- bin.ops[[i]]
	name <- bin.op.strs[[i]]
	for (left.type in type.set) {
		for (right.type in type.set) {
			test.name <- paste("pair-wise vector", left.spec, right.spec, name,
							   left.type, right.type)

			test_that(test.name, {
					  fm.vec1 <- get.vec(left.type, spec.val=left.spec, percent=0.5)
					  fm.vec2 <- get.vec(right.type, spec.val=right.spec, percent=0.5)
					  if (!is.null(fm.vec1) && !is.null(fm.vec2)) {
						  expect_equal(typeof(fm.vec1), left.type)
						  expect_equal(typeof(fm.vec2), right.type)

						  vec1 <- fm.conv.FM2R(fm.vec1)
						  vec2 <- fm.conv.FM2R(fm.vec2)

						  fm.res <- bin.op(fm.vec1, fm.vec2)
						  res <- bin.op(vec1, vec2)

						  expect_equal(res, fm.conv.FM2R(fm.res))
						  expect_equal(typeof(fm.res), typeof(res)) }})
		}
	}
}
}
}

# Test element-wise matrix matrix operation.
for (i in 1:length(bin.ops)) {
	bin.op <- bin.ops[[i]]
	name <- bin.op.strs[[i]]
	for (left.type in type.set) {
		for (right.type in type.set) {
			test.name <- paste("pair-wise matrix", name,
							   left.type, right.type)
			test_that(test.name, {
					  fm.mat1 <- get.mat(left.type)
					  fm.mat2 <- get.mat(right.type)
					  expect_equal(typeof(fm.mat1), left.type)
					  expect_equal(typeof(fm.mat2), right.type)

					  mat1 <- fm.conv.FM2R(fm.mat1)
					  mat2 <- fm.conv.FM2R(fm.mat2)

					  fm.res <- bin.op(fm.mat1, fm.mat2)
					  res <- bin.op(mat1, mat2)

					  expect_equal(res, fm.conv.FM2R(fm.res))
					  expect_equal(typeof(fm.res), typeof(res)) })
		}
	}
}

get.scalar <- function(type)
{
	if (type == "logical")
		return(TRUE)
	else if (type == "integer")
		return(10)
	else if (type == "double")
		return(runif(1))
	else
		stop("a wrong type")
}

# pmin and pmax can take multiple arguments, so we require that all arguments
# have the same class type.
bin.ops1 <- list(`+`, `-`, `*`, `/`, `==`, `!=`, `>`, `>=`, `<`, `<=`, `^`, `|`, `&`,
				 pmin2, pmax2)
bin.op.strs1 <- list("+", "-", "*", "/", "==", "!=", ">", ">=", "<", "<=", "^", "|", "&",
					 "pmin2", "pmax2")

# test element-wise vector scalar operations.
for (i in 1:length(bin.ops1)) {
	bin.op <- bin.ops1[[i]]
	name <- bin.op.strs1[[i]]
	for (left.type in type.set) {
		for (right.type in type.set) {
			test.name <- paste("pair-wise vector element", name,
							   left.type, right.type)

			test_that(test.name, {
					  fm.vec <- get.vec(left.type)
					  expect_equal(typeof(fm.vec), left.type)
					  vec <- fm.conv.FM2R(fm.vec)
					  val <- get.scalar(right.type)

					  fm.res <- bin.op(fm.vec, val)
					  res <- bin.op(vec, val)

					  expect_equal(res, fm.conv.FM2R(fm.res))
					  expect_equal(typeof(fm.res), typeof(res)) })
		}
	}
}

# test element-wise scalar vector operations.
for (i in 1:length(bin.ops1)) {
	bin.op <- bin.ops1[[i]]
	name <- bin.op.strs1[[i]]
	for (left.type in type.set) {
		for (right.type in type.set) {
			test.name <- paste("pair-wise element vector", name,
							   left.type, right.type)

			test_that(test.name, {
					  fm.vec <- get.vec(left.type)
					  expect_equal(typeof(fm.vec), left.type)
					  vec <- fm.conv.FM2R(fm.vec)
					  val <- get.scalar(right.type)

					  fm.res <- bin.op(val, fm.vec)
					  res <- bin.op(val, vec)

					  expect_equal(res, fm.conv.FM2R(fm.res))
					  expect_equal(typeof(fm.res), typeof(res)) })
		}
	}
}

# test element-wise matrix vector operations.
for (i in 1:length(bin.ops1)) {
	bin.op <- bin.ops1[[i]]
	name <- bin.op.strs1[[i]]
	for (left.type in type.set) {
		for (right.type in type.set) {
			test.name <- paste("pair-wise matrix vector", name,
							   left.type, right.type)

			test_that(test.name, {
					  fm.mat <- get.mat(left.type, nrow=20, ncol=100)
					  fm.vec <- get.vec(right.type, len=20)
					  expect_equal(typeof(fm.mat), left.type)
					  expect_equal(typeof(fm.vec), right.type)
					  mat <- fm.conv.FM2R(fm.mat)
					  vec <- fm.conv.FM2R(fm.vec)

					  fm.res <- bin.op(fm.mat, fm.vec)
					  res <- bin.op(mat, vec)

					  expect_equal(res, fm.conv.FM2R(fm.res))
					  expect_equal(typeof(fm.res), typeof(res)) })
		}
	}
}

# test element-wise vector matrix operations.
for (i in 1:length(bin.ops1)) {
	bin.op <- bin.ops1[[i]]
	name <- bin.op.strs1[[i]]
	for (left.type in type.set) {
		for (right.type in type.set) {
			test.name <- paste("pair-wise vector matrix", name,
							   left.type, right.type)

			test_that(test.name, {
					  fm.mat <- get.mat(left.type, nrow=20, ncol=100)
					  fm.vec <- get.vec(right.type, len=20)
					  expect_equal(typeof(fm.mat), left.type)
					  expect_equal(typeof(fm.vec), right.type)
					  mat <- fm.conv.FM2R(fm.mat)
					  vec <- fm.conv.FM2R(fm.vec)

					  fm.res <- bin.op(fm.vec, fm.mat)
					  res <- bin.op(vec, mat)

					  expect_equal(res, fm.conv.FM2R(fm.res))
					  expect_equal(typeof(fm.res), typeof(res)) })
		}
	}
}

# test element-wise matrix Rvector operations.
for (i in 1:length(bin.ops1)) {
	bin.op <- bin.ops1[[i]]
	name <- bin.op.strs1[[i]]
	for (left.type in type.set) {
		for (right.type in type.set) {
			test.name <- paste("pair-wise matrix Rvector", name,
							   left.type, right.type)

			test_that(test.name, {
					  fm.mat <- get.mat(left.type, nrow=20, ncol=100)
					  fm.vec <- get.vec(right.type, len=20)
					  expect_equal(typeof(fm.mat), left.type)
					  expect_equal(typeof(fm.vec), right.type)
					  mat <- fm.conv.FM2R(fm.mat)
					  vec <- fm.conv.FM2R(fm.vec)

					  fm.res <- bin.op(fm.mat, vec)
					  res <- bin.op(mat, vec)

					  expect_equal(res, fm.conv.FM2R(fm.res))
					  expect_equal(typeof(fm.res), typeof(res)) })
		}
	}
}

# test element-wise Rvector matrix operations.
for (i in 1:length(bin.ops1)) {
	bin.op <- bin.ops1[[i]]
	name <- bin.op.strs1[[i]]
	for (left.type in type.set) {
		for (right.type in type.set) {
			test.name <- paste("pair-wise Rvector matrix", name,
							   left.type, right.type)

			test_that(test.name, {
					  fm.mat <- get.mat(left.type, nrow=20, ncol=100)
					  fm.vec <- get.vec(right.type, len=20)
					  expect_equal(typeof(fm.mat), left.type)
					  expect_equal(typeof(fm.vec), right.type)
					  mat <- fm.conv.FM2R(fm.mat)
					  vec <- fm.conv.FM2R(fm.vec)

					  fm.res <- bin.op(vec, fm.mat)
					  res <- bin.op(vec, mat)

					  expect_equal(res, fm.conv.FM2R(fm.res))
					  expect_equal(typeof(fm.res), typeof(res)) })
		}
	}
}

uops <- list(`-`, `!`, abs, sqrt, log, log10, log2, exp, round, ceiling,
			 floor, as.integer, as.numeric)
uop.strs <- list("-", "!", "abs", "sqrt", "log", "log10", "log2", "exp",
				 "round", "ceiling", "floor", "as.integer", "as.numeric")

# Test unary operations.
for (type in type.set) {
for (spec in spec.vals) {
for (i in 1:length(uops)) {
	op <- uops[[i]]
	name <- uop.strs[[i]]
	test_that(paste("test vector", name, type, spec), {
			  fm.vec <- get.vec(type, spec.val=spec, percent=0.5)
			  if (!is.null(fm.vec)) {
				  rvec <- fm.conv.FM2R(fm.vec)
				  fm.res <- op(fm.vec)
				  res <- op(rvec)
				  expect_equal(typeof(fm.res), typeof(res))
				  expect_equal(fm.conv.FM2R(fm.res), res) }})
	test_that(paste("test matrix", name, type, spec), {
			  fm.mat <- get.mat(type, spec.val=spec, percent=0.5)
			  # as.integer and as.numeric in R outputs a vector regardless of
			  # the input object. As to FlashR, these two functions output
			  # an object of the same type as the input object.
			  if (!is.null(fm.mat) && name != "as.integer" && name != "as.numeric") {
				  rmat <- fm.conv.FM2R(fm.mat)
				  fm.res <- op(fm.mat)
				  res <- op(rmat)
				  expect_equal(fm.conv.FM2R(fm.res), res) }})
}
}
}

test_that("test transpose", {
		  mat <- get.mat("double", 20, 100)
		  len <- dim(mat)[1] * dim(mat)[2]
		  fm.t.mat <- t(mat)
		  expect_equal(typeof(fm.t.mat), typeof(mat))
		  mat <- fm.conv.FM2R(mat)
		  t.mat <- t(mat)
		  expect_equal(t.mat, fm.conv.FM2R(fm.t.mat))

		  vec <- get.vec("double")
		  fm.t.mat <- t(vec)
		  expect_true(fm.is.matrix(fm.t.mat))
		  expect_equal(nrow(fm.t.mat), 1)
		  expect_equal(ncol(fm.t.mat), length(vec))
		  expect_equal(as.vector(fm.t.mat), as.vector(vec))
})

# ERROR: boolean indexing fails.
for (type in type.set) {
name <- paste("Get columns from a matrix", type)
test_that(name, {
		  fm.mat <- get.mat(type, nrow=100, ncol=20)
		  mat <- fm.conv.FM2R(fm.mat)
		  fm.sub.mat <- fm.mat[, 3:12]
		  sub.mat <- mat[, 3:12]
		  expect_equal(typeof(fm.sub.mat), typeof(sub.mat))
		  expect_equal(sub.mat, fm.conv.FM2R(fm.sub.mat))

		  fm.mat <- t(fm.mat)
		  mat <- t(mat)
		  fm.sub.mat <- fm.mat[, 3:12]
		  sub.mat <- mat[, 3:12]
		  expect_equal(typeof(fm.sub.mat), typeof(sub.mat))
		  expect_equal(sub.mat, fm.conv.FM2R(fm.sub.mat))
})

name <- paste("Get columns with boolean index from a matrix", type)
test_that(name, {
		  fm.mat <- get.mat(type, 20, 100)
		  mat <- fm.conv.FM2R(fm.mat)
		  fm.index <- fm.runif(ncol(fm.mat)) > 0.5
		  fm.sub.mat <- fm.mat[,fm.index]
		  sub.mat <- mat[,fm.conv.FM2R(fm.index)]
		  expect_equal(typeof(fm.sub.mat), typeof(sub.mat))
		  expect_equal(sub.mat, fm.conv.FM2R(fm.sub.mat))

		  fm.index <- fm.rep.int(FALSE, ncol(fm.mat))
		  fm.sub.mat <- fm.mat[,fm.index]
		  sub.mat <- mat[,fm.conv.FM2R(fm.index)]
		  expect_equal(typeof(fm.sub.mat), typeof(sub.mat))
		  expect_equal(sub.mat, fm.conv.FM2R(fm.sub.mat))

		  fm.index <- fm.runif(ncol(fm.mat) + 1) > 0.5
		  fm.sub.mat <- fm.mat[,fm.index]
		  sub.mat <- mat[,fm.conv.FM2R(fm.index)]
		  expect_equal(typeof(fm.sub.mat), typeof(sub.mat))
		  expect_equal(sub.mat, fm.conv.FM2R(fm.sub.mat))

		  fm.index <- fm.runif(ncol(fm.mat) - 1) > 0.5
		  fm.sub.mat <- fm.mat[,fm.index]
		  sub.mat <- mat[,fm.conv.FM2R(fm.index)]
		  expect_equal(typeof(fm.sub.mat), typeof(sub.mat))
		  expect_equal(sub.mat, fm.conv.FM2R(fm.sub.mat))
})

name <- paste("Get rows from a matrix", type)
test_that(name, {
		  fm.mat <- get.mat(type, nrow=100, ncol=20)
		  mat <- fm.conv.FM2R(fm.mat)
		  fm.sub.mat <- fm.mat[3:12,]
		  sub.mat <- mat[3:12,]
		  expect_equal(typeof(fm.sub.mat), typeof(sub.mat))
		  expect_equal(sub.mat, fm.conv.FM2R(fm.sub.mat))

		  fm.mat <- t(fm.mat)
		  mat <- t(mat)
		  fm.sub.mat <- fm.mat[3:12,]
		  sub.mat <- mat[3:12,]
		  expect_equal(typeof(fm.sub.mat), typeof(sub.mat))
		  expect_equal(sub.mat, fm.conv.FM2R(fm.sub.mat))
})

name <- paste("Get rows with boolean index from a matrix", type)
test_that(name, {
		  fm.mat <- get.mat(type, 20, 100)
		  mat <- fm.conv.FM2R(fm.mat)
		  fm.index <- fm.runif(nrow(fm.mat)) > 0.5
		  fm.sub.mat <- fm.mat[fm.index,]
		  sub.mat <- mat[fm.conv.FM2R(fm.index),]
		  expect_equal(typeof(fm.sub.mat), typeof(sub.mat))
		  expect_equal(sub.mat, fm.conv.FM2R(fm.sub.mat))

		  fm.index <- fm.rep.int(FALSE, nrow(fm.mat))
		  fm.sub.mat <- fm.mat[fm.index,]
		  sub.mat <- mat[fm.conv.FM2R(fm.index),]
		  expect_equal(typeof(fm.sub.mat), typeof(sub.mat))
		  expect_equal(sub.mat, fm.conv.FM2R(fm.sub.mat))

		  fm.index <- fm.runif(nrow(fm.mat) + 1) > 0.5
		  fm.sub.mat <- fm.mat[fm.index,]
		  sub.mat <- mat[fm.conv.FM2R(fm.index),]
		  expect_equal(typeof(fm.sub.mat), typeof(sub.mat))
		  expect_equal(sub.mat, fm.conv.FM2R(fm.sub.mat))

		  fm.index <- fm.runif(nrow(fm.mat) - 1) > 0.5
		  fm.sub.mat <- fm.mat[fm.index,]
		  sub.mat <- mat[fm.conv.FM2R(fm.index),]
		  expect_equal(typeof(fm.sub.mat), typeof(sub.mat))
		  expect_equal(sub.mat, fm.conv.FM2R(fm.sub.mat))
})

name <- paste("Get elements from a matrix", type)
test_that(name, {
		  fm.mat <- get.mat(type, nrow=100, ncol=20)
		  mat <- fm.conv.FM2R(fm.mat)
		  fm.sub.mat <- fm.mat[3:12,3:12]
		  sub.mat <- mat[3:12,3:12]
		  expect_equal(typeof(fm.sub.mat), typeof(sub.mat))
		  expect_equal(sub.mat, fm.conv.FM2R(fm.sub.mat))

		  fm.mat <- t(fm.mat)
		  mat <- t(mat)
		  fm.sub.mat <- fm.mat[3:12,3:12]
		  sub.mat <- mat[3:12,3:12]
		  expect_equal(typeof(fm.sub.mat), typeof(sub.mat))
		  expect_equal(sub.mat, fm.conv.FM2R(fm.sub.mat))
})

name <- paste("Get elements from a vector", type)
test_that(name, {
		  fm.vec <- get.vec(type, 2000)
		  vec <- fm.conv.FM2R(fm.vec)
		  fm.sub.vec <- fm.vec[3:12]
		  sub.vec <- vec[3:12]
		  expect_equal(typeof(fm.sub.vec), typeof(sub.vec))
		  expect_equal(sub.vec, fm.conv.FM2R(fm.sub.vec))
})

name <- paste("Get elements with boolean index from a vector", type)
test_that(name, {
		  fm.vec <- get.vec(type, 2000)
		  vec <- fm.conv.FM2R(fm.vec)
		  fm.index <- fm.runif(length(fm.vec)) > 0.5
		  fm.sub.vec <- fm.vec[fm.index]
		  sub.vec <- vec[fm.conv.FM2R(fm.index)]
		  expect_equal(typeof(fm.sub.vec), typeof(sub.vec))
		  expect_equal(sub.vec, fm.conv.FM2R(fm.sub.vec))

		  fm.index <- fm.rep.int(FALSE, length(fm.vec))
		  fm.sub.vec <- fm.vec[fm.index]
		  sub.vec <- vec[fm.conv.FM2R(fm.index)]
		  expect_equal(typeof(fm.sub.vec), typeof(sub.vec))
		  expect_equal(sub.vec, fm.conv.FM2R(fm.sub.vec))

		  fm.index <- fm.runif(length(fm.vec) + 1) > 0.5
		  fm.sub.vec <- fm.vec[fm.index]
		  sub.vec <- vec[fm.conv.FM2R(fm.index)]
		  expect_equal(typeof(fm.sub.vec), typeof(sub.vec))
		  expect_equal(sub.vec, fm.conv.FM2R(fm.sub.vec))

		  fm.index <- fm.runif(length(fm.vec) - 1) > 0.5
		  fm.sub.vec <- fm.vec[fm.index]
		  sub.vec <- vec[fm.conv.FM2R(fm.index)]
		  expect_equal(typeof(fm.sub.vec), typeof(sub.vec))
		  expect_equal(sub.vec, fm.conv.FM2R(fm.sub.vec))
})

name <- paste("Get head from a matrix", type)
test_that(name, {
		  fm.mat <- get.mat(type, nrow=100, ncol=20)
		  mat <- fm.conv.FM2R(fm.mat)
		  fm.sub.mat <- head(fm.mat)
		  sub.mat <- head(mat)
		  expect_equal(typeof(fm.sub.mat), typeof(sub.mat))
		  expect_equal(sub.mat, fm.conv.FM2R(fm.sub.mat))

		  fm.mat <- t(fm.mat)
		  mat <- t(mat)
		  fm.sub.mat <- head(fm.mat)
		  sub.mat <- head(mat)
		  expect_equal(typeof(fm.sub.mat), typeof(sub.mat))
		  expect_equal(sub.mat, fm.conv.FM2R(fm.sub.mat))
})

name <- paste("Get tail from a matrix", type)
test_that(name, {
		  fm.mat <- get.mat(type, nrow=100, ncol=20)
		  mat <- fm.conv.FM2R(fm.mat)
		  fm.sub.mat <- tail(fm.mat)
		  sub.mat <- tail(mat)
		  expect_equal(typeof(fm.sub.mat), typeof(sub.mat))
		  orig.nrow <- nrow(sub.mat)
		  orig.ncol <- ncol(sub.mat)
		  # tail in R retrieves the last rows from a matrix and still indicates
		  # that the data is retrieved from the tail of the original matrix.
		  # We have to convert it to a vector and back to a matrix to remove
		  # that information.
		  expect_equal(matrix(as.vector(sub.mat), orig.nrow, orig.ncol),
					   fm.conv.FM2R(fm.sub.mat))

		  fm.mat <- t(fm.mat)
		  mat <- t(mat)
		  fm.sub.mat <- tail(fm.mat)
		  sub.mat <- tail(mat)
		  orig.nrow <- nrow(sub.mat)
		  orig.ncol <- ncol(sub.mat)
		  expect_equal(typeof(fm.sub.mat), typeof(sub.mat))
		  expect_equal(matrix(as.vector(sub.mat), orig.nrow, orig.ncol),
					   fm.conv.FM2R(fm.sub.mat))
})

name <- paste("Get head from a vector", type)
test_that(name, {
		  fm.vec <- get.vec(type)
		  vec <- fm.conv.FM2R(fm.vec)
		  fm.sub.vec <- head(fm.vec)
		  sub.vec <- head(vec)
		  expect_equal(typeof(fm.sub.vec), typeof(sub.vec))
		  expect_equal(sub.vec, fm.conv.FM2R(fm.sub.vec))
})

name <- paste("Get tail from a vector", type)
test_that(name, {
		  fm.vec <- get.vec(type)
		  vec <- fm.conv.FM2R(fm.vec)
		  fm.sub.vec <- tail(fm.vec)
		  sub.vec <- tail(vec)
		  expect_equal(typeof(fm.sub.vec), typeof(sub.vec))
		  expect_equal(sub.vec, fm.conv.FM2R(fm.sub.vec))
})
}

agg.ops1 <- list(sum, min, max, any, all)
agg.op.strs1 <- list("sum", "min", "max", "any", "all")
for (na.rm in c(TRUE, FALSE)) {
for (spec1 in spec.vals) {
for (spec2 in spec.vals) {
for (type in type.set) {
	for (i in 1:length(agg.ops1)) {
		agg.op <- agg.ops1[[i]]
		name <- agg.op.strs1[[i]]
		test.name <- paste("aggregate vector", spec1, spec2, name, type, na.rm)

		test_that(test.name, {
				  fm.vec1 <- get.vec(type, spec.val=spec1, percent=0.5)
				  fm.vec2 <- get.vec(type, spec.val=spec2, percent=0.5)
				  if (!is.null(fm.vec1) && !is.null(fm.vec2)) {
					  expect_equal(typeof(fm.vec1), type)
					  expect_equal(typeof(fm.vec2), type)
					  vec1 <- fm.conv.FM2R(fm.vec1)
					  vec2 <- fm.conv.FM2R(fm.vec2)

					  fm.res <- fm.conv.FM2R(agg.op(fm.vec1, na.rm=na.rm))
					  res <- agg.op(vec1, na.rm=na.rm)
					  expect_equal(res, fm.res)
					  expect_equal(typeof(fm.res), typeof(res))

					  fm.res <- fm.conv.FM2R(agg.op(fm.vec1, fm.vec2, na.rm=na.rm))
					  res <- agg.op(vec1, vec2, na.rm=na.rm)
					  expect_equal(res, fm.res)
					  expect_equal(typeof(fm.res), typeof(res)) }})

	}
	test_that(paste("test range", type, spec1, spec2, na.rm), {
			  fm.vec1 <- get.vec(type, spec.val=spec1, percent=0.5)
			  fm.vec2 <- get.vec(type, spec.val=spec2, percent=0.5)
			  if (!is.null(fm.vec1) && !is.null(fm.vec2)) {
				  expect_equal(typeof(fm.vec1), type)
				  expect_equal(typeof(fm.vec2), type)
				  vec1 <- fm.conv.FM2R(fm.vec1)
				  vec2 <- fm.conv.FM2R(fm.vec2)

				  fm.res <- fm.conv.FM2R(range(fm.vec1, na.rm=na.rm))
				  res <- range(vec1, na.rm=na.rm)
				  expect_equal(res, fm.res)
				  expect_equal(typeof(fm.res), typeof(res))

				  fm.res <- fm.conv.FM2R(range(fm.vec1, fm.vec2, na.rm=na.rm))
				  res <- range(vec1, vec2, na.rm=na.rm)
				  expect_equal(res, fm.res)
				  expect_equal(typeof(fm.res), typeof(res)) }})
	test_that(paste("test mean", type, spec1, spec2, na.rm), {
			  fm.vec1 <- get.vec(type, spec.val=spec1, percent=0.5)
			  fm.vec2 <- get.vec(type, spec.val=spec2, percent=0.5)
			  if (!is.null(fm.vec1) && !is.null(fm.vec2)) {
				  expect_equal(typeof(fm.vec1), type)
				  expect_equal(typeof(fm.vec2), type)
				  vec1 <- fm.conv.FM2R(fm.vec1)
				  vec2 <- fm.conv.FM2R(fm.vec2)

				  fm.res <- fm.conv.FM2R(mean(fm.vec1, na.rm=na.rm))
				  res <- mean(vec1, na.rm=na.rm)
				  expect_equal(res, fm.res)
				  expect_equal(typeof(fm.res), typeof(res)) }})

}
}
}
}


for (na.rm in c(TRUE, FALSE)) {
for (i in 1:length(agg.ops1)) {
	agg.op <- agg.ops1[[i]]
	name <- agg.op.strs1[[i]]
	for (type in type.set) {
		test.name <- paste("aggregate matrix", name, type, na.rm)

		test_that(test.name, {
				  fm.mat1 <- get.mat(type)
				  fm.mat2 <- get.mat(type)
				  if (!is.null(fm.mat1) && !is.null(fm.mat2)) {
					  expect_equal(typeof(fm.mat1), type)
					  expect_equal(typeof(fm.mat2), type)
					  mat1 <- fm.conv.FM2R(fm.mat1)
					  mat2 <- fm.conv.FM2R(fm.mat2)

					  fm.res <- fm.conv.FM2R(agg.op(fm.mat1, na.rm=na.rm))
					  res <- agg.op(mat1, na.rm=na.rm)
					  expect_equal(res, fm.res)
					  expect_equal(typeof(fm.res), typeof(res))

					  fm.res <- fm.conv.FM2R(agg.op(fm.mat1, fm.mat2, na.rm=na.rm))
					  res <- agg.op(mat1, mat2, na.rm=na.rm)
					  expect_equal(res, fm.res)
					  expect_equal(typeof(fm.res), typeof(res)) }})
	}
}
}

for (type in type.set) {
for (na.rm in c(TRUE, FALSE)) {
for (spec in spec.vals) {
test_that(paste("pmax on vectors", type, spec, na.rm), {
		  fm.vec1 <- get.vec(type, spec.val=spec, percent=0.5)
		  fm.vec2 <- get.vec(type, spec.val=spec, percent=0.5)
		  fm.vec3 <- get.vec(type, spec.val=spec, percent=0.5)
		  if (!is.null(fm.vec1) && !is.null(fm.vec2) && !is.null(fm.vec3)) {
			  vec1 <- fm.conv.FM2R(fm.vec1)
			  vec2 <- fm.conv.FM2R(fm.vec2)
			  vec3 <- fm.conv.FM2R(fm.vec3)
			  expect_equal(fm.conv.FM2R(pmax(fm.vec1, na.rm=na.rm)), pmax(vec1, na.rm=na.rm))
			  expect_equal(fm.conv.FM2R(pmax(fm.vec1, fm.vec2, na.rm=na.rm)), pmax(vec1, vec2, na.rm=na.rm))
			  expect_equal(fm.conv.FM2R(pmax(fm.vec1, fm.vec2, fm.vec3, na.rm=na.rm)),
						   pmax(vec1, vec2, vec3, na.rm=na.rm))
		  }
})
test_that(paste("pmin on vectors", type, spec, na.rm), {
		  fm.vec1 <- get.vec(type, spec.val=spec, percent=0.5)
		  fm.vec2 <- get.vec(type, spec.val=spec, percent=0.5)
		  fm.vec3 <- get.vec(type, spec.val=spec, percent=0.5)
		  if (!is.null(fm.vec1) && !is.null(fm.vec2) && !is.null(fm.vec3)) {
			  vec1 <- fm.conv.FM2R(fm.vec1)
			  vec2 <- fm.conv.FM2R(fm.vec2)
			  vec3 <- fm.conv.FM2R(fm.vec3)
			  expect_equal(fm.conv.FM2R(pmin(fm.vec1, na.rm=na.rm)), pmin(vec1, na.rm=na.rm))
			  expect_equal(fm.conv.FM2R(pmin(fm.vec1, fm.vec2, na.rm=na.rm)), pmin(vec1, vec2, na.rm=na.rm))
			  expect_equal(fm.conv.FM2R(pmin(fm.vec1, fm.vec2, fm.vec3, na.rm=na.rm)),
						   pmin(vec1, vec2, vec3, na.rm=na.rm))
		  }
})
test_that(paste("pmax/pmin on matrices", type, spec, na.rm), {
		  fm.mat1 <- get.mat(type, spec.val=spec, percent=0.5)
		  fm.mat2 <- get.mat(type, spec.val=spec, percent=0.5)
		  fm.mat3 <- get.mat(type, spec.val=spec, percent=0.5)
		  if (!is.null(fm.mat1) && !is.null(fm.mat2) && !is.null(fm.mat3)) {
			  mat1 <- fm.conv.FM2R(fm.mat1)
			  mat2 <- fm.conv.FM2R(fm.mat2)
			  mat3 <- fm.conv.FM2R(fm.mat3)
			  expect_equal(fm.conv.FM2R(pmax(fm.mat1, na.rm=na.rm)), pmax(mat1, na.rm=na.rm))
			  expect_equal(fm.conv.FM2R(pmax(fm.mat1, fm.mat2, na.rm=na.rm)), pmax(mat1, mat2, na.rm=na.rm))
			  expect_equal(fm.conv.FM2R(pmax(fm.mat1, fm.mat2, fm.mat3, na.rm=na.rm)),
						   pmax(mat1, mat2, mat3, na.rm=na.rm))

			  expect_equal(fm.conv.FM2R(pmin(fm.mat1, na.rm=na.rm)), pmin(mat1, na.rm=na.rm))
			  expect_equal(fm.conv.FM2R(pmin(fm.mat1, fm.mat2, na.rm=na.rm)), pmin(mat1, mat2, na.rm=na.rm))
			  expect_equal(fm.conv.FM2R(pmin(fm.mat1, fm.mat2, fm.mat3, na.rm=na.rm)),
						   pmin(mat1, mat2, mat3, na.rm=na.rm))
		  }
})
}
}
}

for (type in type.set) {
for (na.rm in c(TRUE, FALSE)) {
for (spec in spec.vals) {
test_that(paste("rowSums", type, spec, na.rm), {
		  fm.mat <- get.mat(type, spec.val=spec, percent=0.5)
		  if (!is.null(fm.mat)) {
			  res1 <- fm.conv.FM2R(rowSums(fm.mat, na.rm))
			  res2 <- rowSums(fm.conv.FM2R(fm.mat), na.rm)
			  expect_equal(res1, res2)
		  }
})

test_that(paste("colSums", type, spec, na.rm), {
		  fm.mat <- get.mat(type, spec.val=spec, percent=0.5)
		  if (!is.null(fm.mat)) {
			  res1 <- fm.conv.FM2R(colSums(fm.mat, na.rm))
			  res2 <- colSums(fm.conv.FM2R(fm.mat), na.rm)
			  expect_equal(res1, res2)
		  }
})

test_that(paste("rowMeans", type, spec, na.rm), {
		  fm.mat <- get.mat(type, spec.val=spec, percent=0.5)
		  if (!is.null(fm.mat)) {
			  res1 <- fm.conv.FM2R(rowMeans(fm.mat, na.rm))
			  res2 <- rowMeans(fm.conv.FM2R(fm.mat), na.rm)
			  expect_equal(res1, res2)
		  }
})

test_that(paste("colMeans", type, spec, na.rm), {
		  fm.mat <- get.mat(type, spec.val=spec, percent=0.5)
		  if (!is.null(fm.mat)) {
			  res1 <- fm.conv.FM2R(colMeans(fm.mat, na.rm))
			  res2 <- colMeans(fm.conv.FM2R(fm.mat), na.rm)
			  expect_equal(res1, res2)
		  }
})
}
}
}

test_that("test ifelse", {
		  v1 <- fm.runif(1000)
		  test <- v1 > 0.5
		  res <- ifelse(test, 0, v1)
		  rres <- ifelse(fm.conv.FM2R(test), 0, fm.conv.FM2R(v1))
		  expect_equal(rres, fm.conv.FM2R(res))
		  res <- ifelse(test, fm.rep.int(0, length(v1)), v1)
		  expect_equal(rres, fm.conv.FM2R(res))

		  res <- ifelse(test, v1, 0)
		  rres <- ifelse(fm.conv.FM2R(test), fm.conv.FM2R(v1), 0)
		  expect_equal(rres, fm.conv.FM2R(res))

		  m1 <- fm.runif.matrix(100, 10)
		  test <- m1 > 0.5
		  res <- ifelse(test, 0, m1)
		  rres <- ifelse(fm.conv.FM2R(test), 0, fm.conv.FM2R(m1))
		  expect_equal(rres, fm.conv.FM2R(res))
		  res <- ifelse(test, fm.matrix(0, nrow(m1), ncol(m1)), m1)
		  expect_equal(rres, fm.conv.FM2R(res))

		  res <- ifelse(test, m1, 0)
		  rres <- ifelse(fm.conv.FM2R(test), fm.conv.FM2R(m1), 0)
		  expect_equal(rres, fm.conv.FM2R(res))
})

approx.equal <- function(x1, x2)
{
	fm.conv.FM2R(abs(x1 - x2) < 1e-16)
}

test.with.na <- function(data, fun)
{
	fm.data <- fm.conv.R2FM(data)
	res1 <- fm.conv.FM2R(fun(fm.data))
	res2 <- fm.conv.FM2R(fun(fm.data, na.rm=TRUE))
	if (is.na(res1))
		expect_true(is.na(fun(data)))
	else
		expect_true(fun(data) == res1)
	expect_true(!is.na(res2)
				&& approx.equal(res2, fun(data, na.rm=TRUE)))
}

test_that("Test is.na", {
		  v1 <- runif(1000)
		  v1.int <- as.integer(v1 * 100)
		  v1[v1 < 0.2] <- NA
		  v1.int[v1.int < 20] <- NA
		  expect_equal(is.na(v1), fm.conv.FM2R(is.na(fm.conv.R2FM(v1))))
		  expect_equal(is.na(v1.int), fm.conv.FM2R(is.na(fm.conv.R2FM(v1.int))))

		  m1 <- matrix(v1, 100, 10)
		  m1.int <- matrix(v1.int, 100, 10)
		  expect_equal(is.na(m1), fm.conv.FM2R(is.na(fm.conv.R2FM(m1))))
		  expect_equal(is.na(m1.int), fm.conv.FM2R(is.na(fm.conv.R2FM(m1.int))))

		  v1 <- runif(100) < 0.5
		  v1[1] <- NA
		  test.with.na(v1, all)
		  test.with.na(v1, any)
		  test.with.na(v1, sum)
		  test.with.na(v1, min)
		  test.with.na(v1, max)
		  test.with.na(v1, sd)
		  test.with.na(v1, mean)
		  v1 <- 1:100
		  v1[1] <- NA
		  test.with.na(v1, sum)
		  test.with.na(v1, min)
		  test.with.na(v1, max)
		  test.with.na(v1, sd)
		  test.with.na(v1, mean)
})

test_that("Test is.nan", {
		  v1 <- runif(1000)
		  v1.int <- as.integer(v1 * 100)
		  v1[v1 < 0.2] <- NaN
		  expect_equal(is.nan(v1), fm.conv.FM2R(is.nan(fm.conv.R2FM(v1))))
		  expect_equal(is.nan(v1.int), fm.conv.FM2R(is.nan(fm.conv.R2FM(v1.int))))

		  m1 <- matrix(v1, 100, 10)
		  m1.int <- matrix(v1.int, 100, 10)
		  expect_equal(is.nan(m1), fm.conv.FM2R(is.nan(fm.conv.R2FM(m1))))
		  expect_equal(is.nan(m1.int), fm.conv.FM2R(is.nan(fm.conv.R2FM(m1.int))))

		  v1 <- 1:100
		  v1[1] <- NaN
		  test.with.na(v1, sum)
		  test.with.na(v1, min)
		  test.with.na(v1, max)
		  test.with.na(v1, sd)
		  test.with.na(v1, mean)
})

test_that("Test is.finite", {
		  v1 <- runif(1000)
		  v1.int <- as.integer(v1 * 100)
		  v1[v1 < 0.2] <- Inf
		  expect_equal(is.infinite(v1), fm.conv.FM2R(is.infinite(fm.conv.R2FM(v1))))
		  expect_equal(is.infinite(v1.int), fm.conv.FM2R(is.infinite(fm.conv.R2FM(v1.int))))
		  expect_equal(is.finite(v1), fm.conv.FM2R(is.finite(fm.conv.R2FM(v1))))
		  expect_equal(is.finite(v1.int), fm.conv.FM2R(is.finite(fm.conv.R2FM(v1.int))))

		  m1 <- matrix(v1, 100, 10)
		  m1.int <- matrix(v1.int, 100, 10)
		  expect_equal(is.infinite(m1), fm.conv.FM2R(is.infinite(fm.conv.R2FM(m1))))
		  expect_equal(is.infinite(m1.int), fm.conv.FM2R(is.infinite(fm.conv.R2FM(m1.int))))
		  expect_equal(is.finite(m1), fm.conv.FM2R(is.finite(fm.conv.R2FM(m1))))
		  expect_equal(is.finite(m1.int), fm.conv.FM2R(is.finite(fm.conv.R2FM(m1.int))))
})

test_that("Test matrix info", {
		  fm.mat <- fm.runif.matrix(10, 20)
		  mat <- fm.conv.FM2R(fm.mat)
		  expect_equal(dim(fm.mat), dim(mat))
		  expect_equal(nrow(fm.mat), nrow(mat))
		  expect_equal(ncol(fm.mat), ncol(mat))
		  expect_equal(length(fm.mat), length(mat))
		  expect_equal(typeof(fm.mat), typeof(mat))
		  expect_equal(fm.conv.FM2R(t(fm.mat)), t(mat))

		  fm.vec <- fm.runif(100)
		  vec <- fm.conv.FM2R(fm.vec)
		  expect_equal(dim(fm.vec), dim(vec))
		  expect_equal(nrow(fm.vec), nrow(vec))
		  expect_equal(ncol(fm.vec), ncol(vec))
		  expect_equal(length(fm.vec), length(vec))
		  expect_equal(typeof(fm.vec), typeof(vec))
		  expect_equal(fm.conv.FM2R(t(fm.vec)), t(vec))
})

test_that("test table", {
		  fm.vec <- fm.runif(2000, min=0, max=100)
		  fm.vec <- as.integer(fm.vec)
		  vec <- fm.conv.FM2R(fm.vec)
		  fm.res <- fm.table(fm.vec)
		  res <- as.data.frame(table(vec))
		  expect_equal(as.integer(res$vec), fm.conv.FM2R(fm.res@val + 1))
		  expect_equal(res$Freq, fm.conv.FM2R(fm.res@Freq))
})

cast.type <- function(fm.obj, obj, type)
{
	if (to.type == "double") {
		fm.new.obj <- as.numeric(fm.obj)
		new.obj <- as.numeric(obj)
	}
	else if (to.type == "integer") {
		fm.new.obj <- as.integer(fm.obj)
		new.obj <- as.integer(obj)
	}
	else if (to.type == "logical") {
		fm.new.obj <- as.logical(fm.obj)
		new.obj <- as.logical(obj)
	}
	list(fm.new.obj=fm.new.obj, new.obj=new.obj)
}

for (from.type in c("integer", "double")) {
for (to.type in c("integer", "double")) {
name <- paste("Test casting from", from.type, "to", to.type)
test_that(name, {
		  fm.mat <- get.mat(from.type)
		  mat <- fm.conv.FM2R(fm.mat)
		  res <- cast.type(fm.mat, mat, to.type)
		  fm.new.mat <- res$fm.new.obj
		  new.mat <- matrix(res$new.obj, nrow(mat), ncol(mat))
		  expect_equal(typeof(new.mat), typeof(fm.new.mat))
		  expect_equal(new.mat, fm.conv.FM2R(fm.new.mat))
		  expect_equal(is.numeric(new.mat), is.numeric(fm.new.mat))

		  fm.vec <- get.vec(from.type)
		  vec <- fm.conv.FM2R(fm.vec)
		  res <- cast.type(fm.vec, vec, to.type)
		  fm.new.vec <- res$fm.new.obj
		  new.vec <- res$new.obj
		  expect_equal(typeof(new.vec), typeof(fm.new.vec))
		  expect_equal(new.vec, fm.conv.FM2R(fm.new.vec))
		  expect_equal(is.numeric(new.vec), is.numeric(fm.new.vec))
})
}
}

test_that("test sweep", {
		  m <- fm.runif.matrix(100, 10)
		  v1 <- fm.runif(ncol(m))
		  v2 <- fm.runif(nrow(m))
		  rv1 <- fm.conv.FM2R(v1)
		  rv2 <- fm.conv.FM2R(v2)
		  res <- sweep(m, 1, v2)
		  rres <- sweep(fm.conv.FM2R(m), 1, rv2)
		  expect_equal(fm.conv.FM2R(res), rres)

		  res <- sweep(m, 2, v1)
		  rres <- sweep(fm.conv.FM2R(m), 2, rv1)
		  expect_equal(fm.conv.FM2R(res), rres)
})

test_that("test crossprod", {
		  fm.mat <- fm.runif.matrix(100, 20)
		  fm.mat2 <- fm.runif.matrix(100, 20)
		  mat <- fm.conv.FM2R(fm.mat)
		  mat2 <- fm.conv.FM2R(fm.mat2)

		  fm.res <- crossprod(fm.mat)
		  res <- crossprod(mat)
		  expect_equal(typeof(fm.res), typeof(res))
		  expect_equal(fm.conv.FM2R(fm.res), res)

		  fm.res <- crossprod(fm.mat, fm.mat2)
		  res <- crossprod(mat, mat2)
		  expect_equal(typeof(fm.res), typeof(res))
		  expect_equal(fm.conv.FM2R(fm.res), res)

		  fm.res <- tcrossprod(fm.mat)
		  res <- tcrossprod(mat)
		  expect_equal(typeof(fm.res), typeof(res))
		  expect_equal(fm.conv.FM2R(fm.res), res)

		  fm.res <- tcrossprod(fm.mat, fm.mat2)
		  res <- tcrossprod(mat, mat2)
		  expect_equal(typeof(fm.res), typeof(res))
		  expect_equal(fm.conv.FM2R(fm.res), res)
})

for (type in type.set) {
test_that("which.max", {
		  fm.mat <- get.mat(type, nrow=100, ncol=10)
		  agg.which.max <- fm.create.agg.op(fm.bo.which.max, NULL, "which.max")
		  res1 <- fm.conv.FM2R(fm.agg.mat(fm.mat, 1, agg.which.max))
		  res2 <- apply(fm.conv.FM2R(fm.mat), 1, function(x) which.max(x))
		  expect_equal(res1, res2)
})

test_that("which.min", {
		  fm.mat <- get.mat(type, nrow=100, ncol=10)
		  agg.which.min <- fm.create.agg.op(fm.bo.which.min, NULL, "which.min")
		  res1 <- fm.conv.FM2R(fm.agg.mat(fm.mat, 1, agg.which.min))
		  res2 <- apply(fm.conv.FM2R(fm.mat), 1, function(x) which.min(x))
		  expect_equal(res1, res2)
})
}
