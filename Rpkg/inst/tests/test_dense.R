library(FlashR)

context("Test small matrices")

test_that("Test sequence generator",{
		  fm.vec <- fm.seq.int(1, 10000000, 1)
		  vec <- seq.int(1, 10000000, 1)
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

test_that("convert a FlashMatrixR vector to R vector", {
		  fm.vec <- fm.rep.int(1, 2000)
		  vec <- fm.conv.FM2R(fm.vec)
		  expect_true(is.vector(vec))
		  expect_true(is.atomic(vec))
		  fm.vec <- fm.rep.int(TRUE, 2000)
		  vec <- fm.conv.FM2R(fm.vec)
		  expect_equal(typeof(fm.vec), "logical")
		  expect_true(is.vector(vec))
		  expect_true(is.atomic(vec))
		  expect_equal(vec, fm.conv.FM2R(fm.vec))
})

test_that("convert a R vector to FlashMatrixR vector", {
		  vec <- 1:2000
		  fm.vec <- fm.conv.R2FM(vec)
		  expect_equal(length(fm.vec), length(vec))
		  vec1 <- fm.conv.FM2R(fm.vec)
		  expect_equal(vec, vec1)
		  vec <- runif(1000) < 0.5
		  fm.vec <- fm.conv.R2FM(vec)
		  expect_equal(typeof(fm.vec), "logical")
		  expect_equal(fm.conv.FM2R(fm.vec), vec)
})

test_that("create a column-wise FlashMatrixR matrix", {
		  fm.mat <- fm.matrix(1, 20, 100)
		  expect_equal(dim(fm.mat)[1], 20)
		  expect_equal(dim(fm.mat)[2], 100)
		  mat <- fm.conv.FM2R(fm.mat)
		  expect_equal(nrow(mat), 20)
		  expect_equal(ncol(mat), 100)
		  fm.mat <- fm.matrix(TRUE, 20, 100)
		  expect_equal(typeof(fm.mat), "logical")
})

test_that("convert a R vector to a FlashMatrixR matrix", {
		  mat <- matrix(1:2000, 20, 100)
		  fm.mat <- fm.conv.R2FM(mat)
		  expect_equal(dim(fm.mat)[1], 20)
		  expect_equal(dim(fm.mat)[2], 100)
		  mat1 <- fm.conv.FM2R(fm.mat)
		  expect_equal(mat1, mat)
		  mat <- matrix(runif(2000) < 0.5, 20, 100)
		  fm.mat <- fm.conv.R2FM(mat)
		  expect_equal(typeof(fm.mat), "logical")
})

test.MV1 <- function(fm.mat, fm.vec, mat, vec)
{
	res <- mat %*% vec
	fm.res <- fm.mat %*% fm.vec
	expect_equal(res[, 1], fm.conv.FM2R(fm.res))
	expect_true(fm.is.vector(fm.res))
	expect_equal(length(fm.res), dim(fm.mat)[1])
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

test_that("matrix vector multiply with different type", {
		  mat <- matrix(runif(2000), 20, 100)
		  vec <- runif(100)
		  stopifnot(typeof(mat) == "double")
		  stopifnot(typeof(vec) == "double")
		  test.MV1(fm.conv.R2FM(mat), fm.conv.R2FM(vec), mat, vec)

		  mat <- matrix(1:2000, 20, 100)
		  vec <- 1:100
		  stopifnot(typeof(mat) == "integer")
		  stopifnot(typeof(vec) == "integer")
		  test.MV1(fm.conv.R2FM(mat), fm.conv.R2FM(vec), mat, vec)

		  mat <- matrix(1:2000, 20, 100)
		  vec <- runif(100)
		  stopifnot(typeof(mat) == "integer")
		  stopifnot(typeof(vec) == "double")
		  test.MV1(fm.conv.R2FM(mat), fm.conv.R2FM(vec), mat, vec)

		  mat <- matrix(runif(2000), 20, 100)
		  vec <- 1:100
		  stopifnot(typeof(mat) == "double")
		  stopifnot(typeof(vec) == "integer")
		  test.MV1(fm.conv.R2FM(mat), fm.conv.R2FM(vec), mat, vec)
})

test_that("matrix matrix multiply with different type", {
		  mat1 <- matrix(runif(2000), 20, 100)
		  mat2 <- matrix(runif(2000), 100, 20)
		  stopifnot(typeof(mat1) == "double")
		  stopifnot(typeof(mat2) == "double")
		  test.MM1(fm.conv.R2FM(mat1), fm.conv.R2FM(mat2), mat1, mat2)

		  mat1 <- matrix(1:2000, 20, 100)
		  mat2 <- matrix(1:2000, 100, 20)
		  stopifnot(typeof(mat1) == "integer")
		  stopifnot(typeof(mat2) == "integer")
		  test.MM1(fm.conv.R2FM(mat1), fm.conv.R2FM(mat2), mat1, mat2)

		  mat1 <- matrix(1:2000, 20, 100)
		  mat2 <- matrix(runif(2000), 100, 20)
		  stopifnot(typeof(mat1) == "integer")
		  stopifnot(typeof(mat2) == "double")
		  test.MM1(fm.conv.R2FM(mat1), fm.conv.R2FM(mat2), mat1, mat2)

		  mat1 <- matrix(runif(2000), 20, 100)
		  mat2 <- matrix(1:2000, 100, 20)
		  stopifnot(typeof(mat1) == "double")
		  stopifnot(typeof(mat2) == "integer")
		  test.MM1(fm.conv.R2FM(mat1), fm.conv.R2FM(mat2), mat1, mat2)
})

test_that("column-wise matrix times a vector", {
		  mat <- matrix(1:2000, 20, 100)
		  vec <- 1:ncol(mat)

		  fm.mat <- fm.conv.R2FM(mat, byrow = FALSE)
		  fm.vec <- fm.conv.R2FM(vec)
		  expect_equal(fm.matrix.layout(fm.mat), "col")
		  test.MV1(fm.mat, fm.vec, mat, vec)
})

test_that("Row-wise matrix times a vector", {
		  mat <- matrix(1:2000, 20, 100)
		  vec <- 1:ncol(mat)

		  fm.mat <- fm.conv.R2FM(mat, byrow = TRUE)
		  fm.vec <- fm.conv.R2FM(vec)
		  expect_equal(fm.matrix.layout(fm.mat), "row")
		  test.MV1(fm.mat, fm.vec, mat, vec)
})

test.MM.tmp <- function(left.byrow, right.byrow, left.mat, right.mat)
{
	if (left.byrow)
		left.layout = "row"
	else
		left.layout = "col"
	if (right.byrow)
		right.layout = "row"
	else
		right.layout = "col"
	fm.left <- fm.conv.R2FM(left.mat, byrow = left.byrow)
	expect_equal(fm.matrix.layout(fm.left), left.layout)
	fm.right <- fm.conv.R2FM(right.mat, byrow = right.byrow)
	expect_equal(fm.matrix.layout(fm.right), right.layout)
	test.MM1(fm.left, fm.right, left.mat, right.mat)
}

test_that("matrix multiply: wide vs. tall", {
		  left.mat <- matrix(runif(2000), 20, 100)
		  right.mat <- matrix(runif(2000), 100, 20)
		  # left col-wise, right col-wise
		  test.MM.tmp(FALSE, FALSE, left.mat, right.mat)
		  # left col-wise, right row-wise
		  test.MM.tmp(FALSE, TRUE, left.mat, right.mat)
		  # left row-wise, right col-wise
		  test.MM.tmp(TRUE, FALSE, left.mat, right.mat)
		  # left row-wise, right row-wise
		  test.MM.tmp(TRUE, TRUE, left.mat, right.mat)
})

test_that("matrix multiply: tall vs. small", {
		  left.mat <- matrix(runif(2000), 100, 20)
		  right.mat <- matrix(runif(2000), 20, 10)
		  # left col-wise, right col-wise
		  test.MM.tmp(FALSE, FALSE, left.mat, right.mat)
		  # left col-wise, right row-wise
		  test.MM.tmp(FALSE, TRUE, left.mat, right.mat)
		  # left row-wise, right col-wise
		  test.MM.tmp(TRUE, FALSE, left.mat, right.mat)
		  # left row-wise, right row-wise
		  test.MM.tmp(TRUE, TRUE, left.mat, right.mat)
})

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

type.set <- c("double", "integer", "logical")
bin.ops <- list(`+`, `-`, `*`, `/`, `==`, `!=`, `>`, `>=`, `<`, `<=`,
				pmin, pmax)
bin.op.strs <- list("+", "-", "*", "/", "==", "!=", ">", ">=", "<", "<=",
				"pmin", "pmax")
spec.vals <- list(NULL, "NA", "NaN", "Inf", "-Inf")

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

get.mat <- function(type, nrow=20, ncol=100)
{
	if (type == "integer")
		return(as.integer(fm.seq.matrix(1, 2000, nrow, ncol)))
	else if (type == "double")
		return(fm.runif.matrix(nrow, ncol))
	else if (type == "logical")
		return(fm.seq.matrix(1, 2000, nrow, ncol) < 1000)
	else
		stop("wrong left type")
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
bin.ops1 <- list(`+`, `-`, `*`, `/`, `==`, `!=`, `>`, `>=`, `<`, `<=`,
				 pmin2, pmax2)
bin.op.strs1 <- list("+", "-", "*", "/", "==", "!=", ">", ">=", "<", "<=",
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

uops <- list(abs, sqrt, log, log10, log2, exp, round)
uop.strs <- list("abs", "sqrt", "log", "log10", "log2", "exp", "round")

# Test unary operations on vectors.
for (i in 1:length(uops)) {
	op <- uops[[i]]
	name <- uop.strs[[i]]
	test_that(paste("test vector", name), {
			  fm.vec <- fm.runif(2000)
			  rvec <- fm.conv.FM2R(fm.vec)
			  fm.res <- op(fm.vec)
			  res <- op(rvec)
			  expect_equal(fm.conv.FM2R(fm.res), res) })
}

# Test unary operations on matrices.
for (i in 1:length(uops)) {
	op <- uops[[i]]
	name <- uop.strs[[i]]
	test_that(paste("test matrix", name), {
			  fm.mat <- fm.runif.matrix(20, 100)
			  rmat <- fm.conv.FM2R(fm.mat)
			  fm.res <- op(fm.mat)
			  res <- op(rmat)
			  expect_equal(fm.conv.FM2R(fm.res), res) })
}

test_that("test not", {
		  fm.vec1 <- fm.runif(2000, -1, 1)
		  rvec1 <- fm.conv.FM2R(fm.vec1)
		  fm.vec2 <- fm.runif(2000, -1, 1)
		  rvec2 <- fm.conv.FM2R(fm.vec2)
		  fm.res <- !(fm.vec1 == fm.vec2)
		  res <- !(rvec1 == rvec2)
		  expect_equal(fm.conv.FM2R(fm.res), res)
})

test.t <- function(mat)
{
	len <- dim(mat)[1] * dim(mat)[2]
	fm.t.mat <- t(mat)
	expect_equal(typeof(fm.t.mat), typeof(mat))
	mat <- fm.conv.FM2R(mat)
	t.mat <- t(mat)
	expect_equal(t.mat, fm.conv.FM2R(fm.t.mat))
}

test_that("test transpose", {
		  fm.col.mat <- fm.seq.matrix(1, 2000, 20, 100)
		  test.t(fm.col.mat)
		  fm.col.mat <- fm.runif.matrix(20, 100)
		  test.t(fm.col.mat)
})

test.get.col <- function(fm.mat)
{
	mat <- fm.conv.FM2R(fm.mat)
	fm.sub.mat <- fm.mat[, 3:12]
	sub.mat <- mat[, 3:12]
	expect_equal(sub.mat, fm.conv.FM2R(fm.sub.mat))
}

test_that("Get columns from a column-wise matrix", {
		  fm.mat <- fm.seq.matrix(1, 2000, 100, 20)
		  test.get.col(fm.mat)
		  fm.mat <- fm.runif.matrix(20, 100) < 0.5
		  test.get.col(fm.mat)
})

test.get.row <- function(fm.mat)
{
	mat <- fm.conv.FM2R(fm.mat)
	fm.sub.mat <- fm.mat[3:12,]
	sub.mat <- mat[3:12,]
	expect_equal(sub.mat, fm.conv.FM2R(fm.sub.mat))
}

test_that("Get rows from a column-wise matrix", {
		  fm.mat <- fm.seq.matrix(1, 2000, 100, 20)
		  test.get.row(fm.mat)
		  fm.mat <- fm.runif.matrix(20, 100) < 0.5
		  test.get.row(fm.mat)
})

test_that("read/write a dense matrix", {
		  fm.mat <- fm.runif.matrix(100, 20)
		  mat <- fm.conv.FM2R(fm.mat)
		  expect_true(fm.write.obj(fm.mat, "test.mat"))
		  fm.mat1 <- fm.read.obj("test.mat")
		  expect_equal(fm.matrix.layout(fm.mat1), "col")
		  expect_equal(mat, fm.conv.FM2R(fm.mat1))
})

test_that("read/write a transposed dense matrix", {
		  fm.mat <- fm.runif.matrix(100, 20)
		  expect_true(fm.write.obj(t(fm.mat), "test.mat"))
		  fm.mat1 <- fm.read.obj("test.mat")
		  expect_equal(fm.matrix.layout(fm.mat1), "row")
		  expect_equal(fm.conv.FM2R(t(fm.mat)), fm.conv.FM2R(fm.mat1))
})

test_that("read/write a dense submatrix", {
		  fm.mat <- fm.runif.matrix(100, 20)
		  expect_true(fm.write.obj(fm.get.cols(fm.mat, 3:12), "test.mat"))
		  fm.mat1 <- fm.read.obj("test.mat")
		  expect_equal(fm.matrix.layout(fm.mat1), "col")
		  expect_equal(fm.conv.FM2R(fm.mat[,3:12]), fm.conv.FM2R(fm.mat1))
})

test_that("read/write a transposed dense submatrix", {
		  fm.mat <- fm.runif.matrix(100, 20)
		  expect_true(fm.write.obj(t(fm.get.cols(fm.mat, 3:12)), "test.mat"))
		  fm.mat1 <- fm.read.obj("test.mat")
		  expect_equal(fm.matrix.layout(fm.mat1), "row")
		  expect_equal(fm.conv.FM2R(t(fm.mat[,3:12])), fm.conv.FM2R(fm.mat1))
})

agg.ops1 <- list(sum, min, max)
agg.op.strs1 <- list("sum", "min", "max")
for (na.rm in c(TRUE, FALSE)) {
for (spec1 in spec.vals) {
for (spec2 in spec.vals) {
for (i in 1:length(agg.ops1)) {
	agg.op <- agg.ops1[[i]]
	name <- agg.op.strs1[[i]]
	for (type in type.set) {
		test.name <- paste("aggregate vector", spec1, spec2, name, type, na.rm)

		test_that(test.name, {
				  fm.vec1 <- get.vec(type, spec.val=spec1, percent=0.5)
				  fm.vec2 <- get.vec(type, spec.val=spec2, percent=0.5)
				  if (!is.null(fm.vec1) && !is.null(fm.vec2)) {
					  expect_equal(typeof(fm.vec1), type)
					  expect_equal(typeof(fm.vec2), type)
					  vec1 <- fm.conv.FM2R(fm.vec1)
					  vec2 <- fm.conv.FM2R(fm.vec2)

					  fm.res <- agg.op(fm.vec1, na.rm=na.rm)
					  res <- agg.op(vec1, na.rm=na.rm)
					  expect_equal(res, fm.res)
					  expect_equal(typeof(fm.res), typeof(res))

					  fm.res <- agg.op(fm.vec1, fm.vec2, na.rm=na.rm)
					  res <- agg.op(vec1, vec2, na.rm=na.rm)
					  expect_equal(res, fm.res)
					  expect_equal(typeof(fm.res), typeof(res)) }})
	}
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

					  fm.res <- agg.op(fm.mat1, na.rm=na.rm)
					  res <- agg.op(mat1, na.rm=na.rm)
					  expect_equal(res, fm.res)
					  expect_equal(typeof(fm.res), typeof(res))

					  fm.res <- agg.op(fm.mat1, fm.mat2, na.rm=na.rm)
					  res <- agg.op(mat1, mat2, na.rm=na.rm)
					  expect_equal(res, fm.res)
					  expect_equal(typeof(fm.res), typeof(res)) }})
	}
}
}

test_that("pmax on vectors", {
		  fm.vec1 <- fm.runif(2000)
		  fm.vec2 <- fm.runif(2000)
		  fm.vec3 <- fm.runif(2000)
		  vec1 <- fm.conv.FM2R(fm.vec1)
		  vec2 <- fm.conv.FM2R(fm.vec2)
		  vec3 <- fm.conv.FM2R(fm.vec3)
		  expect_equal(fm.conv.FM2R(pmax(fm.vec1)), pmax(vec1))
		  expect_equal(fm.conv.FM2R(pmax(fm.vec1, fm.vec2)), pmax(vec1, vec2))
		  expect_equal(fm.conv.FM2R(pmax(fm.vec1, fm.vec2, fm.vec3)),
					   pmax(vec1, vec2, vec3))
})

test_that("pmin on vectors", {
		  fm.vec1 <- fm.runif(2000)
		  fm.vec2 <- fm.runif(2000)
		  fm.vec3 <- fm.runif(2000)
		  vec1 <- fm.conv.FM2R(fm.vec1)
		  vec2 <- fm.conv.FM2R(fm.vec2)
		  vec3 <- fm.conv.FM2R(fm.vec3)
		  expect_equal(fm.conv.FM2R(pmin(fm.vec1)), pmin(vec1))
		  expect_equal(fm.conv.FM2R(pmin(fm.vec1, fm.vec2)), pmin(vec1, vec2))
		  expect_equal(fm.conv.FM2R(pmin(fm.vec1, fm.vec2, fm.vec3)),
					   pmin(vec1, vec2, vec3))
})

test_that("pmax/pmin on matrices", {
		  fm.mat1 <- fm.runif.matrix(20, 100)
		  mat1 <- fm.conv.FM2R(fm.mat1)
		  fm.mat2 <- fm.runif.matrix(20, 100)
		  mat2 <- fm.conv.FM2R(fm.mat2)
		  fm.mat3 <- fm.runif.matrix(20, 100)
		  mat3 <- fm.conv.FM2R(fm.mat3)
		  expect_equal(fm.conv.FM2R(pmax(fm.mat1)), pmax(mat1))
		  expect_equal(fm.conv.FM2R(pmax(fm.mat1, fm.mat2)), pmax(mat1, mat2))
		  expect_equal(fm.conv.FM2R(pmax(fm.mat1, fm.mat2, fm.mat3)),
					   pmax(mat1, mat2, mat3))

		  expect_equal(fm.conv.FM2R(pmin(fm.mat1)), pmin(mat1))
		  expect_equal(fm.conv.FM2R(pmin(fm.mat1, fm.mat2)), pmin(mat1, mat2))
		  expect_equal(fm.conv.FM2R(pmin(fm.mat1, fm.mat2, fm.mat3)),
					   pmin(mat1, mat2, mat3))
})


test_that("rowSums", {
		  fm.mat <- fm.runif.matrix(100, 10)
		  res1 <- fm.conv.FM2R(rowSums(fm.mat))
		  res2 <- rowSums(fm.conv.FM2R(fm.mat))
		  expect_equal(res1, res2)
})

test_that("colSums", {
		  fm.mat <- fm.runif.matrix(100, 10)
		  res1 <- fm.conv.FM2R(colSums(fm.mat))
		  res2 <- colSums(fm.conv.FM2R(fm.mat))
		  expect_equal(res1, res2)
})

test_that("rowMeans", {
		  fm.mat <- fm.runif.matrix(100, 10)
		  res1 <- fm.conv.FM2R(rowMeans(fm.mat))
		  res2 <- rowMeans(fm.conv.FM2R(fm.mat))
		  expect_equal(res1, res2)
})

test_that("colMeans", {
		  fm.mat <- fm.runif.matrix(100, 10)
		  res1 <- fm.conv.FM2R(colMeans(fm.mat))
		  res2 <- colMeans(fm.conv.FM2R(fm.mat))
		  expect_equal(res1, res2)
})

test_that("which.max", {
		  fm.mat <- fm.runif.matrix(100, 10)
		  agg.which.max <- fm.create.agg.op(fm.bo.which.max, NULL, "which.max")
		  res1 <- fm.conv.FM2R(fm.agg.mat(fm.mat, 1, agg.which.max))
		  res2 <- apply(fm.conv.FM2R(fm.mat), 1, function(x) which.max(x))
		  expect_equal(res1, res2)
})

test_that("which.min", {
		  fm.mat <- fm.runif.matrix(100, 10)
		  agg.which.min <- fm.create.agg.op(fm.bo.which.min, NULL, "which.min")
		  res1 <- fm.conv.FM2R(fm.agg.mat(fm.mat, 1, agg.which.min))
		  res2 <- apply(fm.conv.FM2R(fm.mat), 1, function(x) which.min(x))
		  expect_equal(res1, res2)
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
		  rres[1,] <- colSums(rmat[rlabels == 0,])
		  rres[2,] <- colSums(rmat[rlabels == 1,])
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
		  rres[,1] <- rowSums(rmat[,rlabels == 0])
		  rres[,2] <- rowSums(rmat[,rlabels == 1])
		  expect_equal(fm.conv.FM2R(res), rres)
})

test_that("test sweep", {
		  m <- fm.runif.matrix(100, 10)
		  v1 <- fm.runif(ncol(m))
		  v2 <- fm.runif(nrow(m))
		  rv1 <- fm.conv.FM2R(v1)
		  rv2 <- fm.conv.FM2R(v2)
		  res <- sweep(m, 1, v2)
		  rres <- sweep(fm.conv.FM2R(m), 1, v2)
		  expect_equal(fm.conv.FM2R(res), rres)

		  res <- sweep(m, 2, v1)
		  rres <- sweep(fm.conv.FM2R(m), 2, v1)
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

test_that("test ifelse", {
		  v1 <- fm.runif(1000)
		  test <- v1 > 0.5
		  res <- ifelse(test, 0, v1)
		  rres <- ifelse(fm.conv.FM2R(test), 0, fm.conv.FM2R(v1))
		  expect_equal(rres, fm.conv.FM2R(res))

		  res <- ifelse(test, v1, 0)
		  rres <- ifelse(fm.conv.FM2R(test), fm.conv.FM2R(v1), 0)
		  expect_equal(rres, fm.conv.FM2R(res))

		  m1 <- fm.runif.matrix(100, 10)
		  test <- m1 > 0.5
		  res <- ifelse(test, 0, m1)
		  rres <- ifelse(fm.conv.FM2R(test), 0, fm.conv.FM2R(m1))
		  expect_equal(rres, fm.conv.FM2R(res))

		  res <- ifelse(test, m1, 0)
		  rres <- ifelse(fm.conv.FM2R(test), fm.conv.FM2R(m1), 0)
		  expect_equal(rres, fm.conv.FM2R(res))
})

approx.equal <- function(x1, x2)
{
	abs(x1 - x2) < 1e-16
}

test.with.na <- function(data, fun)
{
	fm.data <- fm.conv.R2FM(data)
	res1 <- fun(fm.data)
	res2 <- fun(fm.data, na.rm=TRUE)
	expect_true(is.na(res1))
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
