library(FlashGraphR)
#fg.set.conf("flash-graph/conf/run_test.txt")

print("Test sequence generator")
fm.vec <- fm.seq.int(1, 10000000, 1)
vec <- seq.int(1, 10000000, 1)
stopifnot(sum(fm.conv.FM2R(fm.vec) == vec) == fm.length(fm.vec))
stopifnot(fm.typeof(fm.vec) == typeof(vec))

print("Test random generator")
fm.vec <- fm.runif(10000000, -1, 1)
vec <- fm.conv.FM2R(fm.vec)
stopifnot(sum(vec >= -1) == fm.length(fm.vec))
stopifnot(sum(vec <= 1) == fm.length(fm.vec))
stopifnot(fm.typeof(fm.vec) == typeof(vec))

# Test a FlashMatrixR vector.
print("convert a FlashMatrixR vector to R vector")
fm.vec <- fm.rep.int(1, 2000)
vec <- fm.conv.FM2R(fm.vec)
stopifnot(is.vector(vec))
stopifnot(is.atomic(vec))

print("convert a R vector to FlashMatrixR vector")
vec <- 1:2000
fm.vec <- fm.conv.R2FM(vec)
stopifnot(fm.length(fm.vec) == length(vec))
vec1 <- fm.conv.FM2R(fm.vec)
stopifnot(sum(vec == vec1) == length(vec))

# Test a FlashMatrixR matrix
print("create a column-wise FlashMatrixR matrix")
fm.mat <- fm.matrix(fm.rep.int(1, 2000), 20, 100)
stopifnot(fm.nrow(fm.mat) == 20)
stopifnot(fm.ncol(fm.mat) == 100)
mat <- fm.conv.FM2R(fm.mat)
stopifnot(nrow(mat) == 20)
stopifnot(ncol(mat) == 100)

print("convert a R vector to a FlashMatrixR matrix")
mat <- matrix(1:2000, 20, 100)
fm.mat <- fm.conv.R2FM(mat)
stopifnot(fm.nrow(fm.mat) == 20)
stopifnot(fm.ncol(fm.mat) == 100)
mat1 <- fm.conv.FM2R(fm.mat)
stopifnot(sum(mat1 == mat) == fm.ncol(fm.mat) * fm.nrow(fm.mat))

test.MV1 <- function(fm.mat, fm.vec, res)
{
	fm.res <- fm.multiply(fm.mat, fm.vec)
	stopifnot(sum(res[, 1] == fm.conv.FM2R(fm.res)) == fm.length(fm.res))
	stopifnot(fm.is.vector(fm.res))
	stopifnot(fm.length(fm.res) == fm.nrow(fm.mat))
	# TODO
#	stopifnot(fm.typeof(fm.res) == typeof(res))
}

# Test a dense matrix times a vector
test.MV <- function(mat, vec)
{
	res <- mat %*% vec
	cat("matrix:", typeof(mat), ", vector:", typeof(vec), ", res:", typeof(res), "\n")
	# TODO
	#print("column-wise matrix times a vector")
	#fm.mat <- fm.conv.R2FM(mat, byrow = FALSE)
	#fm.vec <- fm.conv.R2FM(vec)
	#stopifnot(fm.matrix.layout(fm.mat) == "col")
	#test.MV1(fm.mat, fm.vec, res)

	print("Row-wise matrix times a vector")
	fm.mat <- fm.conv.R2FM(mat, byrow = TRUE)
	fm.vec <- fm.conv.R2FM(vec)
	stopifnot(fm.matrix.layout(fm.mat) == "row")
	test.MV1(fm.mat, fm.vec, res)
}

mat <- matrix(1:2000, 20, 100)
vec <- 1:ncol(mat)
test.MV(mat, vec)

mat <- matrix(1:2000, 20, 100)
vec <- seq.int(1, ncol(mat), 1)
test.MV(mat, vec)

mat <- matrix(seq.int(1, 2000, 1), 20, 100)
vec <- 1:ncol(mat)
test.MV(mat, vec)

mat <- matrix(seq.int(1, 2000, 1), 20, 100)
vec <- seq.int(1, ncol(mat), 1)
test.MV(mat, vec)

test.MM1 <- function(fm.left, fm.right, res)
{
	fm.res <- fm.multiply(fm.left, fm.right)
	stopifnot(sum(res == fm.conv.FM2R(fm.res)) == fm.nrow(fm.res) * fm.ncol(fm.res))
	stopifnot(fm.nrow(fm.res) == fm.nrow(fm.left))
	stopifnot(fm.ncol(fm.res) == fm.ncol(fm.right))
	# TODO
#	stopifnot(fm.typeof(fm.res) == typeof(res))
}

test.MM <- function(left.type, right.type)
{
	if (left.type == "integer")
		left.mat <- matrix(1:2000, 20, 100)
	else
		left.mat <- matrix(seq.int(1, 2000, 1), 20, 100)
	if (right.type == "integer")
		right.mat <- matrix(1:2000, 100, 20)
	else
		right.mat <- matrix(seq.int(1, 2000, 1), 100, 20)

	print("a wide dense matrix times a tall dense matrix")
	res <- left.mat %*% right.mat
	cat("left:", typeof(left.mat), ", right:", typeof(right.mat), ", res:", typeof(res), "\n")

	test.MM.tmp <- function(left.byrow, right.byrow)
	{
		cat("left:", nrow(left.mat), ncol(left.mat), ", right:",
			nrow(right.mat), ncol(right.mat), "\n")
		if (left.byrow)
			left.layout = "row"
		else
			left.layout = "col"
		if (right.byrow)
			right.layout = "row"
		else
			right.layout = "col"
		fm.left <- fm.conv.R2FM(left.mat, byrow = left.byrow)
		stopifnot(fm.matrix.layout(fm.left) == left.layout)
		fm.right <- fm.conv.R2FM(right.mat, byrow = right.byrow)
		stopifnot(fm.matrix.layout(fm.right) == right.layout)
		test.MM1(fm.left, fm.right, res)
	}

	#print("Column-wise matrix times a column-wise matrix")
	#test.MM.tmp(FALSE, FALSE)

	#print("Column-wise matrix times a row-wise matrix")
	#test.MM.tmp(FALSE, TRUE)

	print("Row-wise matrix times a column-wise matrix")
	test.MM.tmp(TRUE, FALSE)

	#print("Row-wise matrix times a row-wise matrix")
	#test.MM.tmp(TRUE, TRUE)

	print("a tall dense matrix times a wide dense matrix")
	if (left.type == "integer")
		left.mat <- matrix(1:2000, 100, 20)
	else
		left.mat <- matrix(seq.int(1, 2000, 1), 100, 20)
	if (right.type == "integer")
		right.mat <- matrix(1:2000, 20, 10)
	else
		right.mat <- matrix(seq.int(1, 2000, 1), 20, 10)
	res <- left.mat %*% right.mat

	print("Column-wise matrix times a column-wise matrix")
	test.MM.tmp(FALSE, FALSE)

	print("Column-wise matrix times a row-wise matrix")
	test.MM.tmp(FALSE, TRUE)

	print("Row-wise matrix times a column-wise matrix")
	test.MM.tmp(TRUE, FALSE)

	#print("Row-wise matrix times a row-wise matrix")
	#test.MM.tmp(TRUE, TRUE)
}

test.MM("integer", "integer")
test.MM("integer", "double")
test.MM("double", "integer")
test.MM("double", "double")

test.basic.op1 <- function(fm.o1, fm.o2, ro1, ro2)
{
	fm.res <- fm.o1 + fm.o2
	res <- ro1 + ro2
	cat("add: left:", typeof(ro1), ", right:", typeof(ro2), ", res:", typeof(res), "\n")
	stopifnot(sum(res == fm.conv.FM2R(fm.res)) == length(res))
	stopifnot(fm.typeof(fm.res) == typeof(res))

	fm.res <- fm.o1 - fm.o2
	res <- ro1 - ro2
	cat("sub left:", typeof(ro1), ", right:", typeof(ro2), ", res:", typeof(res), "\n")
	stopifnot(sum(res == fm.conv.FM2R(fm.res)) == length(res))
	stopifnot(fm.typeof(fm.res) == typeof(res))

	fm.res <- fm.o1 * fm.o2
	res <- ro1 * ro2
	cat("mul left:", typeof(ro1), ", right:", typeof(ro2), ", res:", typeof(res), "\n")
	stopifnot(sum(res == fm.conv.FM2R(fm.res)) == length(res))
	stopifnot(fm.typeof(fm.res) == typeof(res))

	fm.res <- fm.o1 / fm.o2
	res <- ro1 / ro2
	cat("div left:", typeof(ro1), ", right:", typeof(ro2), ", res:", typeof(res), "\n")
	stopifnot(sum(res == fm.conv.FM2R(fm.res)) == length(res))
	stopifnot(fm.typeof(fm.res) == typeof(res))
}

test.basic.op <- function(fm.o1, fm.o2)
{
	ro1 <- fm.conv.FM2R(fm.o1)
	ro2 <- fm.conv.FM2R(fm.o2)

	fm.res <- fm.o1 + fm.o2
	res <- ro1 + ro2
	cat("add left:", typeof(ro1), ", right:", typeof(ro2), ", res:", typeof(res), "\n")
	stopifnot(sum(res == fm.conv.FM2R(fm.res)) == length(res))
	stopifnot(fm.typeof(fm.res) == typeof(res))

	fm.res <- fm.o1 - fm.o2
	res <- ro1 - ro2
	cat("sub left:", typeof(ro1), ", right:", typeof(ro2), ", res:", typeof(res), "\n")
	stopifnot(sum(res == fm.conv.FM2R(fm.res)) == length(res))
	stopifnot(fm.typeof(fm.res) == typeof(res))

	fm.res <- fm.o1 * fm.o2
	res <- ro1 * ro2
	cat("mul left:", typeof(ro1), ", right:", typeof(ro2), ", res:", typeof(res), "\n")
	stopifnot(sum(res == fm.conv.FM2R(fm.res)) == length(res))
	stopifnot(fm.typeof(fm.res) == typeof(res))

	fm.res <- fm.o1 / fm.o2
	res <- ro1 / ro2
	cat("div left:", typeof(ro1), ", right:", typeof(ro2), ", res:", typeof(res), "\n")
	stopifnot(sum(res == fm.conv.FM2R(fm.res)) == length(res))
	stopifnot(fm.typeof(fm.res) == typeof(res))

	fm.res <- fm.pmin2(fm.o1, fm.o2)
	res <- pmin(ro1, ro2)
	cat("min left:", typeof(ro1), ", right:", typeof(ro2), ", res:", typeof(res), "\n")
	stopifnot(sum(res == fm.conv.FM2R(fm.res)) == length(res))
	stopifnot(fm.typeof(fm.res) == typeof(res))

	fm.res <- fm.pmax2(fm.o1, fm.o2)
	res <- pmax(ro1, ro2)
	cat("max left:", typeof(ro1), ", right:", typeof(ro2), ", res:", typeof(res), "\n")
	stopifnot(sum(res == fm.conv.FM2R(fm.res)) == length(res))
	stopifnot(fm.typeof(fm.res) == typeof(res))
}

print("pair-wise integer-integer vector")
fm.vec1 <- fm.conv.R2FM(1:2000)
fm.vec2 <- fm.conv.R2FM(1:2000)
test.basic.op(fm.vec1, fm.vec2)

print("pair-wise double-integer vector")
fm.vec1 <- fm.seq.int(1, 2000, 1)
fm.vec2 <- fm.conv.R2FM(1:2000)
test.basic.op(fm.vec1, fm.vec2)

print("pair-wise integer-double vector")
fm.vec1 <- fm.conv.R2FM(1:2000)
fm.vec2 <- fm.seq.int(1, 2000, 1)
test.basic.op(fm.vec1, fm.vec2)

print("pair-wise double-double vector")
fm.vec1 <- fm.seq.int(1, 2000, 1)
fm.vec2 <- fm.seq.int(1, 2000, 1)
test.basic.op(fm.vec1, fm.vec2)

print("pair-wise integer-integer matrix")
fm.mat1 <- fm.matrix(fm.conv.R2FM(1:2000), 20, 100)
fm.mat2 <- fm.matrix(fm.conv.R2FM(1:2000), 20, 100)
test.basic.op(fm.mat1, fm.mat2)

print("pair-wise double-integer matrix")
fm.mat1 <- fm.matrix(fm.seq.int(1, 2000, 1), 20, 100)
fm.mat2 <- fm.matrix(fm.conv.R2FM(1:2000), 20, 100)
test.basic.op(fm.mat1, fm.mat2)

print("pair-wise integer-double matrix")
fm.mat1 <- fm.matrix(fm.conv.R2FM(1:2000), 20, 100)
fm.mat2 <- fm.matrix(fm.seq.int(1, 2000, 1), 20, 100)
test.basic.op(fm.mat1, fm.mat2)

print("pair-wise double-double matrix")
fm.mat1 <- fm.matrix(fm.seq.int(1, 2000, 1), 20, 100)
fm.mat2 <- fm.matrix(fm.seq.int(1, 2000, 1), 20, 100)
test.basic.op(fm.mat1, fm.mat2)

int.v <- (1:3)[1]

fm.vec <- fm.conv.R2FM(1:2000)
rvec <- fm.conv.FM2R(fm.vec)
print("vector-element integer-integer")
test.basic.op1(fm.vec, int.v, rvec, int.v)
print("vector-element integer-double")
test.basic.op1(fm.vec, 1.0, rvec, 1.0)

fm.vec <- fm.seq.int(1, 2000, 1)
rvec <- fm.conv.FM2R(fm.vec)
print("vector-element double-int")
test.basic.op1(fm.vec, int.v, rvec, int.v)
print("vector-element double-double")
test.basic.op1(fm.vec, 1.0, rvec, 1.0)

fm.vec <- fm.conv.R2FM(1:2000)
rvec <- fm.conv.FM2R(fm.vec)
print("element-vector integer-integer")
test.basic.op1(int.v, fm.vec, int.v, rvec)
print("element-vector integer-double")
test.basic.op1(1.0, fm.vec, 1.0, rvec)

fm.vec <- fm.seq.int(1, 2000, 1)
rvec <- fm.conv.FM2R(fm.vec)
print("element-vector double-int")
test.basic.op1(int.v, fm.vec, int.v, rvec)
print("element-vector double-double")
test.basic.op1(1.0, fm.vec, 1.0, rvec)

print("Transpose a matrix")
fm.mat <- fm.matrix(fm.seq.int(1, 2000, 1), 20, 100)
fm.t.mat <- fm.t(fm.mat)
stopifnot(fm.matrix.layout(fm.t.mat) == "row")
mat <- fm.conv.FM2R(fm.mat)
t.mat <- t(mat)
stopifnot(sum(t.mat == fm.conv.FM2R(fm.t.mat)) == fm.nrow(fm.mat) * fm.ncol(fm.mat))
print("basic operations on a transposed matrix")
test.basic.op(fm.t.mat, fm.matrix(fm.seq.int(1, 2000, 1), 100, 20, TRUE))
print("multiplication on a transposed matrix")
fm.vec <- fm.seq.int(1, 20, 1)
test.MV1(fm.t.mat, fm.vec, fm.conv.FM2R(fm.t.mat) %*% fm.conv.FM2R(fm.vec))
fm.right.mat <- fm.matrix(fm.seq.int(1, 200, 1), 20, 10)
test.MM1(fm.t.mat, fm.right.mat, fm.conv.FM2R(fm.t.mat) %*% fm.conv.FM2R(fm.right.mat))

# TODO test transpose a sparse matrix.

# Test on a sparse matrix
print("load sparse matrix")
fg <- fg.load.graph("data/wiki.adj-v4", "data/wiki.index-v4")
fm <- fm.get.matrix(fg)
stopifnot(!fm.is.sym(fm))
stopifnot(fm.is.sparse(fm))

print("SpMV")
fm.vec <- fm.rep.int(1, fm.ncol(fm))
fm.res <- fm.multiply(fm, fm.vec)
stopifnot(fm.sum(fm.res) == fg.ecount(fg))

# Test SpMM
# TODO

# Test on a very large dense matrix.
# convert a large dense matrix in R to FlashMatrixR.
