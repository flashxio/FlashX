library(FlashR)
#fg.set.conf("flash-graph/conf/run_test.txt")

print("Test sequence generator")
fm.vec <- fm.seq.int(1, 10000000, 1)
vec <- seq.int(1, 10000000, 1)
stopifnot(sum(fm.conv.FM2R(fm.vec) == vec) == length(fm.vec))
stopifnot(typeof(fm.vec) == typeof(vec))

print("Test random generator")
fm.vec <- fm.runif(10000000, -1, 1)
vec <- fm.conv.FM2R(fm.vec)
stopifnot(sum(vec >= -1) == length(fm.vec))
stopifnot(sum(vec <= 1) == length(fm.vec))
stopifnot(typeof(fm.vec) == typeof(vec))

# Test a FlashMatrixR vector.
print("convert a FlashMatrixR vector to R vector")
fm.vec <- fm.rep.int(1, 2000)
vec <- fm.conv.FM2R(fm.vec)
stopifnot(is.vector(vec))
stopifnot(is.atomic(vec))
fm.vec <- fm.rep.int(TRUE, 2000)
vec <- fm.conv.FM2R(fm.vec)
stopifnot(typeof(fm.vec) == "logical")
stopifnot(is.vector(vec))
stopifnot(is.atomic(vec))
stopifnot(sum(vec == fm.vec) == length(vec))

print("convert a R vector to FlashMatrixR vector")
vec <- 1:2000
fm.vec <- fm.conv.R2FM(vec)
stopifnot(length(fm.vec) == length(vec))
vec1 <- fm.conv.FM2R(fm.vec)
stopifnot(sum(vec == vec1) == length(vec))
vec <- runif(1000) < 0.5
fm.vec <- fm.conv.R2FM(vec)
stopifnot(typeof(fm.vec) == "logical")
stopifnot(sum(fm.vec == vec) == length(vec))

# Test a FlashMatrixR matrix
print("create a column-wise FlashMatrixR matrix")
fm.mat <- fm.matrix(fm.rep.int(1, 2000), 20, 100)
stopifnot(dim(fm.mat)[1] == 20)
stopifnot(dim(fm.mat)[2] == 100)
mat <- fm.conv.FM2R(fm.mat)
stopifnot(nrow(mat) == 20)
stopifnot(ncol(mat) == 100)
fm.mat <- fm.matrix(fm.rep.int(TRUE, 2000), 20, 100)
stopifnot(typeof(fm.mat) == "logical")

print("convert a R vector to a FlashMatrixR matrix")
mat <- matrix(1:2000, 20, 100)
fm.mat <- fm.conv.R2FM(mat)
stopifnot(dim(fm.mat)[1] == 20)
stopifnot(dim(fm.mat)[2] == 100)
mat1 <- fm.conv.FM2R(fm.mat)
stopifnot(sum(mat1 == mat) == dim(fm.mat)[2] * dim(fm.mat)[1])
mat <- matrix(runif(2000) < 0.5, 20, 100)
fm.mat <- fm.conv.R2FM(mat)
stopifnot(typeof(fm.mat) == "logical")

test.MV1 <- function(fm.mat, fm.vec, res)
{
	fm.res <- fm.multiply(fm.mat, fm.vec)
	stopifnot(sum(res[, 1] == fm.res) == length(fm.res))
	stopifnot(fm.is.vector(fm.res))
	stopifnot(length(fm.res) == dim(fm.mat)[1])
	stopifnot(typeof(fm.res) == typeof(res))
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
	stopifnot(sum(res == fm.res) == dim(fm.res)[1] * dim(fm.res)[2])
	stopifnot(dim(fm.res)[1] == dim(fm.left)[1])
	stopifnot(dim(fm.res)[2] == dim(fm.right)[2])
	stopifnot(typeof(fm.res) == typeof(res))
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
	stopifnot(sum(res == fm.res) == length(res))
	stopifnot(typeof(fm.res) == typeof(res))

	fm.res <- fm.o1 - fm.o2
	res <- ro1 - ro2
	cat("sub left:", typeof(ro1), ", right:", typeof(ro2), ", res:", typeof(res), "\n")
	stopifnot(sum(res == fm.res) == length(res))
	stopifnot(typeof(fm.res) == typeof(res))

	fm.res <- fm.o1 * fm.o2
	res <- ro1 * ro2
	cat("mul left:", typeof(ro1), ", right:", typeof(ro2), ", res:", typeof(res), "\n")
	stopifnot(sum(res == fm.res) == length(res))
	stopifnot(typeof(fm.res) == typeof(res))

	fm.res <- fm.o1 / fm.o2
	res <- ro1 / ro2
	cat("div left:", typeof(ro1), ", right:", typeof(ro2), ", res:", typeof(res), "\n")
	stopifnot(sum(res == fm.res) == length(res))
	stopifnot(typeof(fm.res) == typeof(res))

	fm.res <- fm.o1 == fm.o2
	res <- ro1 == ro2
	cat("eq left:", typeof(ro1), ", right:", typeof(ro2), ", res:", typeof(res), "\n")
	stopifnot(sum(res == fm.res) == length(res))
	stopifnot(typeof(fm.res) == typeof(res))

	fm.res <- fm.o1 != fm.o2
	res <- ro1 != ro2
	cat("neq left:", typeof(ro1), ", right:", typeof(ro2), ", res:", typeof(res), "\n")
	stopifnot(sum(res == fm.res) == length(res))
	stopifnot(typeof(fm.res) == typeof(res))

	fm.res <- fm.o1 > fm.o2
	res <- ro1 > ro2
	cat("gt left:", typeof(ro1), ", right:", typeof(ro2), ", res:", typeof(res), "\n")
	stopifnot(sum(res == fm.res) == length(res))
	stopifnot(typeof(fm.res) == typeof(res))

	fm.res <- fm.o1 >= fm.o2
	res <- ro1 >= ro2
	cat("ge left:", typeof(ro1), ", right:", typeof(ro2), ", res:", typeof(res), "\n")
	stopifnot(sum(res == fm.res) == length(res))
	stopifnot(typeof(fm.res) == typeof(res))

	fm.res <- fm.o1 < fm.o2
	res <- ro1 < ro2
	cat("lt left:", typeof(ro1), ", right:", typeof(ro2), ", res:", typeof(res), "\n")
	stopifnot(sum(res == fm.res) == length(res))
	stopifnot(typeof(fm.res) == typeof(res))

	fm.res <- fm.o1 <= fm.o2
	res <- ro1 <= ro2
	cat("le left:", typeof(ro1), ", right:", typeof(ro2), ", res:", typeof(res), "\n")
	stopifnot(sum(res == fm.res) == length(res))
	stopifnot(typeof(fm.res) == typeof(res))
}

test.basic.op <- function(fm.o1, fm.o2)
{
	ro1 <- fm.conv.FM2R(fm.o1)
	ro2 <- fm.conv.FM2R(fm.o2)

	fm.res <- fm.o1 + fm.o2
	res <- ro1 + ro2
	cat("add left:", typeof(ro1), ", right:", typeof(ro2), ", res:", typeof(res), "\n")
	stopifnot(sum(res == fm.res) == length(res))
	stopifnot(typeof(fm.res) == typeof(res))

	fm.res <- fm.o1 - fm.o2
	res <- ro1 - ro2
	cat("sub left:", typeof(ro1), ", right:", typeof(ro2), ", res:", typeof(res), "\n")
	stopifnot(sum(res == fm.res) == length(res))
	stopifnot(typeof(fm.res) == typeof(res))

	fm.res <- fm.o1 * fm.o2
	res <- ro1 * ro2
	cat("mul left:", typeof(ro1), ", right:", typeof(ro2), ", res:", typeof(res), "\n")
	stopifnot(sum(res == fm.res) == length(res))
	stopifnot(typeof(fm.res) == typeof(res))

	fm.res <- fm.o1 / fm.o2
	res <- ro1 / ro2
	cat("div left:", typeof(ro1), ", right:", typeof(ro2), ", res:", typeof(res), "\n")
	stopifnot(sum(res == fm.res) == length(res))
	stopifnot(typeof(fm.res) == typeof(res))

	fm.res <- fm.o1 == fm.o2
	res <- ro1 == ro2
	cat("eq left:", typeof(ro1), ", right:", typeof(ro2), ", res:", typeof(res), "\n")
	stopifnot(sum(res == fm.res) == length(res))
	stopifnot(typeof(fm.res) == typeof(res))

	fm.res <- fm.o1 != fm.o2
	res <- ro1 != ro2
	cat("neq left:", typeof(ro1), ", right:", typeof(ro2), ", res:", typeof(res), "\n")
	stopifnot(sum(res == fm.res) == length(res))
	stopifnot(typeof(fm.res) == typeof(res))

	fm.res <- fm.o1 > fm.o2
	res <- ro1 > ro2
	cat("gt left:", typeof(ro1), ", right:", typeof(ro2), ", res:", typeof(res), "\n")
	stopifnot(sum(res == fm.res) == length(res))
	stopifnot(typeof(fm.res) == typeof(res))

	fm.res <- fm.o1 >= fm.o2
	res <- ro1 >= ro2
	cat("ge left:", typeof(ro1), ", right:", typeof(ro2), ", res:", typeof(res), "\n")
	stopifnot(sum(res == fm.res) == length(res))
	stopifnot(typeof(fm.res) == typeof(res))

	fm.res <- fm.o1 < fm.o2
	res <- ro1 < ro2
	cat("lt left:", typeof(ro1), ", right:", typeof(ro2), ", res:", typeof(res), "\n")
	stopifnot(sum(res == fm.res) == length(res))
	stopifnot(typeof(fm.res) == typeof(res))

	fm.res <- fm.o1 <= fm.o2
	res <- ro1 <= ro2
	cat("le left:", typeof(ro1), ", right:", typeof(ro2), ", res:", typeof(res), "\n")
	stopifnot(sum(res == fm.res) == length(res))
	stopifnot(typeof(fm.res) == typeof(res))

	fm.res <- pmin(fm.o1, fm.o2)
	res <- pmin(ro1, ro2)
	cat("min left:", typeof(ro1), ", right:", typeof(ro2), ", res:", typeof(res), "\n")
	stopifnot(sum(res == fm.res) == length(res))
	stopifnot(typeof(fm.res) == typeof(res))

	fm.res <- pmax(fm.o1, fm.o2)
	res <- pmax(ro1, ro2)
	cat("max left:", typeof(ro1), ", right:", typeof(ro2), ", res:", typeof(res), "\n")
	stopifnot(sum(res == fm.res) == length(res))
	stopifnot(typeof(fm.res) == typeof(res))
}

print("pair-wise integer-integer vector")
fm.vec1 <- fm.conv.R2FM(1:2000)
fm.vec2 <- fm.conv.R2FM(1:2000)
test.basic.op(fm.vec1, fm.vec2)
print("pair-wise integer-boolean vector")
test.basic.op(fm.vec1, fm.vec2 < 1000)
print("pair-wise boolean-integer vector")
test.basic.op(fm.vec1 < 1000, fm.vec2)

print("pair-wise double-integer vector")
fm.vec1 <- fm.seq.int(1, 2000, 1)
fm.vec2 <- fm.conv.R2FM(1:2000)
test.basic.op(fm.vec1, fm.vec2)
print("pair-wise double-boolean vector")
test.basic.op(fm.vec1, fm.vec2 < 1000)

print("pair-wise integer-double vector")
fm.vec1 <- fm.conv.R2FM(1:2000)
fm.vec2 <- fm.seq.int(1, 2000, 1)
test.basic.op(fm.vec1, fm.vec2)
print("pair-wise boolean-double vector")
test.basic.op(fm.vec1 < 1000, fm.vec2)

print("pair-wise double-double vector")
fm.vec1 <- fm.seq.int(1, 2000, 1)
fm.vec2 <- fm.seq.int(1, 2000, 1)
test.basic.op(fm.vec1, fm.vec2)

print("pair-wise double-double random vector")
fm.vec1 <- fm.runif(2000)
fm.vec2 <- fm.runif(2000)
test.basic.op(fm.vec1, fm.vec2)

print("pair-wise integer-integer matrix")
fm.mat1 <- fm.matrix(fm.conv.R2FM(1:2000), 20, 100)
fm.mat2 <- fm.matrix(fm.conv.R2FM(1:2000), 20, 100)
test.basic.op(fm.mat1, fm.mat2)
print("pair-wise integer-boolean matrix")
test.basic.op(fm.mat1, fm.mat2 < 1000)
print("pair-wise boolean-integer matrix")
test.basic.op(fm.mat1 < 1000, fm.mat2)

print("pair-wise double-integer matrix")
fm.mat1 <- fm.matrix(fm.seq.int(1, 2000, 1), 20, 100)
fm.mat2 <- fm.matrix(fm.conv.R2FM(1:2000), 20, 100)
test.basic.op(fm.mat1, fm.mat2)
print("pair-wise double-boolean matrix")
test.basic.op(fm.mat1, fm.mat2 < 1000)

print("pair-wise integer-double matrix")
fm.mat1 <- fm.matrix(fm.conv.R2FM(1:2000), 20, 100)
fm.mat2 <- fm.matrix(fm.seq.int(1, 2000, 1), 20, 100)
test.basic.op(fm.mat1, fm.mat2)
print("pair-wise boolean-double matrix")
test.basic.op(fm.mat1 < 1000,  fm.mat2)

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
print("vector-element integer-boolean")
test.basic.op1(fm.vec, TRUE, rvec, TRUE)

fm.vec <- fm.seq.int(1, 2000, 1)
rvec <- fm.conv.FM2R(fm.vec)
print("vector-element double-int")
test.basic.op1(fm.vec, int.v, rvec, int.v)
print("vector-element double-double")
test.basic.op1(fm.vec, 1.0, rvec, 1.0)
print("vector-element double-boolean")
test.basic.op1(fm.vec, TRUE, rvec, TRUE)

fm.vec <- fm.conv.R2FM(1:2000)
rvec <- fm.conv.FM2R(fm.vec)
print("element-vector integer-integer")
test.basic.op1(int.v, fm.vec, int.v, rvec)
print("element-vector integer-double")
test.basic.op1(1.0, fm.vec, 1.0, rvec)
print("element-vector boolean-double")
test.basic.op1(FALSE, fm.vec, FALSE, rvec)

fm.vec <- fm.seq.int(1, 2000, 1)
rvec <- fm.conv.FM2R(fm.vec)
print("element-vector double-int")
test.basic.op1(int.v, fm.vec, int.v, rvec)
print("element-vector double-double")
test.basic.op1(1.0, fm.vec, 1.0, rvec)
print("element-vector boolean-double")
test.basic.op1(FALSE, fm.vec, FALSE, rvec)

test.all.forms <- function(FUN)
{
	fm1 <- fm.matrix(fm.runif(1000), 100, 10)
	fm2 <- fm.matrix(fm.runif(1000), 100, 10)
	fm.res <- FUN(fm1, fm2)
	stopifnot(sum(fm.conv.FM2R(FUN(fm1, fm2))
				  != FUN(fm.conv.FM2R(fm1), fm.conv.FM2R(fm2))) == 0)
	fmv1 <- fm.runif(1000)
	fmv2 <- fm.runif(1000)
	stopifnot(sum(fm.conv.FM2R(FUN(fmv1, fmv2))
				  != FUN(fm.conv.FM2R(fmv1), fm.conv.FM2R(fmv2))) == 0)
	vcol <- fm.runif(nrow(fm1))
	stopifnot(sum(fm.conv.FM2R(FUN(fm1, vcol))
				  != FUN(fm.conv.FM2R(fm1), fm.conv.FM2R(vcol))) == 0)
	stopifnot(sum(fm.conv.FM2R(FUN(vcol, fm1))
				  != FUN(fm.conv.FM2R(vcol), fm.conv.FM2R(fm1))) == 0)
	stopifnot(sum(fm.conv.FM2R(FUN(fm1, fm.conv.FM2R(fm2)))
				  != FUN(fm.conv.FM2R(fm1), fm.conv.FM2R(fm2))) == 0)
	stopifnot(sum(fm.conv.FM2R(FUN(fm.conv.FM2R(fm2), fm1))
				  != FUN(fm.conv.FM2R(fm2), fm.conv.FM2R(fm1))) == 0)
	stopifnot(sum(fm.conv.FM2R(FUN(fm1, 1))
				  != FUN(fm.conv.FM2R(fm1), 1)) == 0)
	stopifnot(sum(fm.conv.FM2R(FUN(1, fm1))
				  != FUN(1, fm.conv.FM2R(fm1))) == 0)
	stopifnot(sum(fm.conv.FM2R(FUN(fm1, fm.conv.FM2R(vcol)))
				  != FUN(fm.conv.FM2R(fm1), fm.conv.FM2R(vcol))) == 0)
	stopifnot(sum(fm.conv.FM2R(FUN(fm.conv.FM2R(vcol), fm1))
				  != FUN(fm.conv.FM2R(vcol), fm.conv.FM2R(fm1))) == 0)
	stopifnot(sum(fm.conv.FM2R(FUN(fmv1, 1))
				  != FUN(fm.conv.FM2R(fmv1), 1)) == 0)
	stopifnot(sum(fm.conv.FM2R(FUN(1, fmv1))
				  != FUN(1, fm.conv.FM2R(fmv1))) == 0)
	stopifnot(sum(fm.conv.FM2R(FUN(fmv1, fm.conv.FM2R(fmv2)))
				  != FUN(fm.conv.FM2R(fmv1), fm.conv.FM2R(fmv2))) == 0)
	stopifnot(sum(fm.conv.FM2R(FUN(fm.conv.FM2R(fmv2), fmv1))
				  != FUN(fm.conv.FM2R(fmv2), fm.conv.FM2R(fmv1))) == 0)
}

test.all.forms(function(x, y) x + y)
test.all.forms(function(x, y) x - y)
test.all.forms(function(x, y) x * y)
test.all.forms(function(x, y) x == y)
test.all.forms(function(x, y) x != y)
test.all.forms(function(x, y) x > y)
test.all.forms(function(x, y) x >= y)
test.all.forms(function(x, y) x <= y)
test.all.forms(function(x, y) x < y)
test.all.forms(function(x, y) pmin2(x, y))
test.all.forms(function(x, y) pmax2(x, y))
#test.all.forms(function(x, y) x / y)

fm.vec1 <- fm.runif(2000, -1, 1)
rvec1 <- fm.conv.FM2R(fm.vec1)
print("abs")
fm.res <- abs(fm.vec1)
res <- abs(rvec1)
stopifnot(sum(fm.res == res) == length(res))
fm.res <- sqrt(fm.res + 1)
res <- sqrt(res + 1)
stopifnot(sum(fm.res == res) == length(res))
fm.vec2 <- fm.runif(2000, -1, 1)
rvec2 <- fm.conv.FM2R(fm.vec2)
stopifnot(sum((!(fm.vec1 == fm.vec2)) == !(rvec1 == rvec2)) == length(rvec1))

test.matrix <- function(fm.mat)
{
	if (fm.matrix.layout(fm.mat) == "row")
		byrow <- TRUE
	else
		byrow <- FALSE
	m.nrow <- dim(fm.mat)[1]
	m.ncol <- dim(fm.mat)[2]
	len <- m.nrow * m.ncol
	mat <- fm.conv.FM2R(fm.mat)
	print("basic operations on a matrix")
	test.basic.op(fm.mat, fm.matrix(fm.seq.int(1, len, 1), m.nrow, m.ncol, byrow))
	print("multiplication on a matrix")
	fm.vec <- fm.seq.int(1, m.ncol, 1)
	test.MV1(fm.mat, fm.vec, fm.conv.FM2R(fm.mat) %*% fm.conv.FM2R(fm.vec))
	fm.right.mat <- fm.matrix(fm.seq.int(1, m.ncol * 10, 1), m.ncol, 10)
	test.MM1(fm.mat, fm.right.mat, fm.conv.FM2R(fm.mat) %*% fm.conv.FM2R(fm.right.mat))
}

test.t.matrix <- function(fm.col.mat)
{
	stopifnot(fm.matrix.layout(fm.col.mat) == "col")
	len <- dim(fm.col.mat)[1] * dim(fm.col.mat)[2]
	fm.t.mat <- t(fm.col.mat)
	stopifnot(typeof(fm.t.mat) == typeof(fm.col.mat))
	t.nrow <- dim(fm.t.mat)[1]
	t.ncol <- dim(fm.t.mat)[2]
	stopifnot(fm.matrix.layout(fm.t.mat) == "row")
	mat <- fm.conv.FM2R(fm.col.mat)
	t.mat <- t(mat)
	stopifnot(sum(t.mat == fm.conv.FM2R(fm.t.mat)) == len)
	print("basic operations on a transposed matrix")
	test.basic.op(fm.t.mat, fm.matrix(fm.seq.int(1, len, 1), t.nrow, t.ncol, TRUE))
	print("multiplication on a transposed matrix")
	fm.vec <- fm.seq.int(1, t.ncol, 1)
	test.MV1(fm.t.mat, fm.vec, fm.conv.FM2R(fm.t.mat) %*% fm.conv.FM2R(fm.vec))
	fm.right.mat <- fm.matrix(fm.seq.int(1, t.ncol * 10, 1), t.ncol, 10)
	test.MM1(fm.t.mat, fm.right.mat, fm.conv.FM2R(fm.t.mat) %*% fm.conv.FM2R(fm.right.mat))
}

print("Transpose a matrix")
fm.mat <- fm.matrix(fm.seq.int(1, 2000, 1), 20, 100)
test.t.matrix(fm.mat)
fm.mat <- fm.matrix(fm.runif(2000) < 0.5, 20, 100)
test.t.matrix(fm.mat)

#print("Set columns of a column-wise matrix")
#fm.mat <- fm.matrix(fm.seq.int(1, 2000, 1), 100, 20)
#mat <- fm.conv.FM2R(fm.mat)
#fm.mat2 <- fm.matrix(fm.seq.int(1001, 1200, 1), 100, 2)
#mat2 <- fm.conv.FM2R(fm.mat2)
#fm.set.cols(fm.mat, 1:2, fm.mat2)
#mat[,1:2] <- mat2
#stopifnot(sum(mat == fm.conv.FM2R(fm.mat)) == length(mat))

print("Get a submatrix from a column-wise matrix")
fm.mat <- fm.matrix(fm.seq.int(1, 2000, 1), 100, 20)
mat <- fm.conv.FM2R(fm.mat)
fm.sub.mat <- fm.get.cols(fm.mat, 3:12)
sub.mat <- mat[, 3:12]
stopifnot(sum(sub.mat == fm.conv.FM2R(fm.sub.mat)) == length(sub.mat))
test.matrix(fm.sub.mat)
test.t.matrix(fm.sub.mat)
fm.mat <- fm.matrix(fm.runif(2000) < 0.5, 20, 100)
fm.sub.mat <- fm.get.cols(fm.mat, 3:12)
stopifnot(typeof(fm.sub.mat) == typeof(fm.mat))

print("read/write a dense matrix")
fm.mat <- fm.matrix(fm.seq.int(1, 2000, 1), 100, 20)
mat <- fm.conv.FM2R(fm.mat)
stopifnot(fm.write.obj(fm.mat, "test.mat"))
fm.mat1 <- fm.read.obj("test.mat")
stopifnot(fm.matrix.layout(fm.mat1) == "col")
stopifnot(sum(mat == fm.conv.FM2R(fm.mat1)) == length(mat))

print("read/write a transposed dense matrix")
stopifnot(fm.write.obj(t(fm.mat), "test.mat"))
fm.mat1 <- fm.read.obj("test.mat")
stopifnot(fm.matrix.layout(fm.mat1) == "row")
stopifnot(sum(t(mat) == fm.conv.FM2R(fm.mat1)) == length(mat))

print("read/write a dense submatrix")
stopifnot(fm.write.obj(fm.get.cols(fm.mat, 3:12), "test.mat"))
fm.mat1 <- fm.read.obj("test.mat")
stopifnot(fm.matrix.layout(fm.mat1) == "col")
stopifnot(sum(mat[,3:12] == fm.conv.FM2R(fm.mat1)) == length(mat[,3:12]))

print("read/write a transposed dense submatrix")
stopifnot(fm.write.obj(t(fm.get.cols(fm.mat, 3:12)), "test.mat"))
fm.mat1 <- fm.read.obj("test.mat")
stopifnot(fm.matrix.layout(fm.mat1) == "row")
stopifnot(sum(t(mat[,3:12]) == fm.conv.FM2R(fm.mat1)) == length(mat[,3:12]))

print("sum on a vector")
fm.vec1 <- fm.runif(2000)
fm.vec2 <- fm.runif(2000)
vec1 <- fm.conv.FM2R(fm.vec1)
vec2 <- fm.conv.FM2R(fm.vec2)
stopifnot(sum(fm.vec1) == sum(vec1))
stopifnot(sum(fm.vec1, fm.vec2) == sum(vec1, vec2))

print("min on a vector")
stopifnot(min(fm.vec1) == min(vec1))
stopifnot(min(fm.vec1, fm.vec2) == min(vec1, vec2))
vec <- 1:2000
fm.vec <- fm.conv.R2FM(vec)
stopifnot(min(fm.vec) == min(vec))

print("max on a vector")
stopifnot(max(fm.vec1) == max(vec1))
stopifnot(max(fm.vec1, fm.vec2) == max(vec1, vec2))
vec <- 1:2000
fm.vec <- fm.conv.R2FM(vec)
stopifnot(max(fm.vec) == max(vec))

print("pmax on vectors")
fm.vec3 <- fm.runif(2000)
vec3 <- fm.conv.FM2R(fm.vec3)
stopifnot(sum(fm.conv.FM2R(pmax(fm.vec1)) == pmax(vec1)) == length(vec1))
stopifnot(sum(fm.conv.FM2R(pmax(fm.vec1, fm.vec2)) == pmax(vec1, vec2)) == length(vec1))
stopifnot(sum(fm.conv.FM2R(pmax(fm.vec1, fm.vec2, fm.vec3))
			  == pmax(vec1, vec2, vec3)) == length(vec1))

print("pmin on vectors")
stopifnot(sum(fm.conv.FM2R(pmin(fm.vec1)) == pmin(vec1)) == length(vec1))
stopifnot(sum(fm.conv.FM2R(pmin(fm.vec1, fm.vec2)) == pmin(vec1, vec2)) == length(vec1))
stopifnot(sum(fm.conv.FM2R(pmin(fm.vec1, fm.vec2, fm.vec3))
			  == pmin(vec1, vec2, vec3)) == length(vec1))

print("sum on a matrix")
fm.mat <- fm.matrix(fm.seq.int(1, 2000, 1), 20, 100)
rmat <- fm.conv.FM2R(fm.mat)
stopifnot(sum(fm.mat) == sum(rmat))
stopifnot(sum(t(fm.mat)) == sum(rmat))

print("pmax on matrices")
fm.mat1 <- fm.matrix(fm.runif(2000), 20, 100)
mat1 <- fm.conv.FM2R(fm.mat1)
fm.mat2 <- fm.matrix(fm.runif(2000), 20, 100)
mat2 <- fm.conv.FM2R(fm.mat2)
fm.mat3 <- fm.matrix(fm.runif(2000), 20, 100)
mat3 <- fm.conv.FM2R(fm.mat3)
stopifnot(sum(fm.conv.FM2R(pmax(fm.mat1)) == pmax(mat1)) == length(mat1))
stopifnot(sum(fm.conv.FM2R(pmax(fm.mat1, fm.mat2)) == pmax(mat1, mat2)) == length(mat1))
stopifnot(sum(fm.conv.FM2R(pmax(fm.mat1, fm.mat2, fm.mat3))
			  == pmax(mat1, mat2, mat3)) == length(mat1))

print("pmin on matrices")
stopifnot(sum(fm.conv.FM2R(pmin(fm.mat1)) == pmin(mat1)) == length(mat1))
stopifnot(sum(fm.conv.FM2R(pmin(fm.mat1, fm.mat2)) == pmin(mat1, mat2)) == length(mat1))
stopifnot(sum(fm.conv.FM2R(pmin(fm.mat1, fm.mat2, fm.mat3))
			  == pmin(mat1, mat2, mat3)) == length(mat1))

print("rowSums")
fm.mat <- fm.matrix(fm.runif(1000), 100, 10)
res1 <- fm.conv.FM2R(rowSums(fm.mat))
res2 <- rowSums(fm.conv.FM2R(fm.mat))
stopifnot(sum(res1 == res2) == length(res1))
print("colSums")
res1 <- fm.conv.FM2R(colSums(fm.mat))
res2 <- colSums(fm.conv.FM2R(fm.mat))
stopifnot(sum(res1 == res2) == length(res1))
print("rowMeans")
res1 <- fm.conv.FM2R(rowMeans(fm.mat))
res2 <- rowMeans(fm.conv.FM2R(fm.mat))
stopifnot(sum(res1 == res2) == length(res1))
print("colMeans")
res1 <- fm.conv.FM2R(colMeans(fm.mat))
res2 <- colMeans(fm.conv.FM2R(fm.mat))
stopifnot(sum(res1 == res2) == length(res1))
print("log")
res1 <- fm.conv.FM2R(log(fm.mat))
res2 <- log(fm.conv.FM2R(fm.mat))
stopifnot(sum(res1 == res2) == length(res1))
print("log10")
res1 <- fm.conv.FM2R(log10(fm.mat))
res2 <- log10(fm.conv.FM2R(fm.mat))
stopifnot(sum(res1 == res2) == length(res1))
print("log2")
res1 <- fm.conv.FM2R(log2(fm.mat))
res2 <- log2(fm.conv.FM2R(fm.mat))
stopifnot(sum(res1 == res2) == length(res1))
print("exp")
res1 <- fm.conv.FM2R(exp(fm.mat))
res2 <- exp(fm.conv.FM2R(fm.mat))
stopifnot(sum(abs(res1 - res2) < 1e-15) == length(res1))
print("round")
res1 <- fm.conv.FM2R(round(fm.mat))
res2 <- round(fm.conv.FM2R(fm.mat))
stopifnot(sum(res1 == res2) == length(res1))

print("which.max")
agg.which.max <- fm.create.agg.op(fm.bo.which.max, NULL, "which.max")
res1 <- fm.conv.FM2R(fm.agg.mat(fm.mat, 1, agg.which.max))
res2 <- apply(fm.conv.FM2R(fm.mat), 1, function(x) which.max(x))
stopifnot(sum(res1 == res2) == length(res1))
print("which.min")
agg.which.min <- fm.create.agg.op(fm.bo.which.min, NULL, "which.min")
res1 <- fm.conv.FM2R(fm.agg.mat(fm.mat, 1, agg.which.min))
res2 <- apply(fm.conv.FM2R(fm.mat), 1, function(x) which.min(x))
stopifnot(sum(res1 == res2) == length(res1))

print("test groupby")
m <- fm.matrix(fm.runif(1000), 100, 10)
v <- floor(fm.runif(100))
labels <- fm.as.factor(fm.as.integer(v), max(v) + 1)
agg.sum <- fm.create.agg.op(fm.bo.add, fm.bo.add, "sum")
res <- fm.groupby(m, 2, labels, agg.sum)
res2 <- fm.agg.mat(m, 2, agg.sum)
stopifnot(sum(fm.as.vector(res) == res2) == length(res2))

v <- floor(fm.runif(100, min=0, max=2))
labels <- fm.as.factor(fm.as.integer(v), max(v) + 1)
agg.sum <- fm.create.agg.op(fm.bo.add, fm.bo.add, "sum")
res <- fm.groupby(m, 2, labels, agg.sum)
rmat <- fm.conv.FM2R(m)
rlabels <- fm.conv.FM2R(v)
rres <- matrix(nrow=2, ncol=10)
rres[1,] <- colSums(rmat[rlabels == 0,])
rres[2,] <- colSums(rmat[rlabels == 1,])
stopifnot(sum(fm.conv.FM2R(res) == rres) == length(rres))

print("test getting rows")
fm.mat <- fm.matrix(fm.seq.int(1, 2000, 1), 100, 20)
rmat <- fm.conv.FM2R(fm.mat)
stopifnot(sum(fm.conv.FM2R(fm.get.rows(fm.mat, 1:10)) == rmat[1:10,])
		  == 10*ncol(rmat))
stopifnot(sum(fm.conv.FM2R(fm.get.cols(fm.mat, 2:11)) == rmat[,2:11])
		  == 10*nrow(rmat))

m <- fm.matrix(fm.runif(1000), 100, 10)
v1 <- fm.runif(ncol(m))
v2 <- fm.runif(nrow(m))
rv1 <- fm.conv.FM2R(v1)
rv2 <- fm.conv.FM2R(v2)
print("test mapply rows")
res <- fm.mapply.row(m, v1, fm.bo.add)
rres <- apply(fm.conv.FM2R(m), 1, function(x) x + rv1)
stopifnot(sum(fm.conv.FM2R(res) == t(rres)) == length(rres))
print("test mapply cols")
res <- fm.mapply.col(m, v2, fm.bo.add)
rres <- apply(fm.conv.FM2R(m), 2, function(x) x + rv2)
stopifnot(sum(fm.conv.FM2R(res) == rres) == length(rres))

m <- fm.matrix(fm.runif(1000), 100, 10)
m2 <- fm.matrix(fm.runif(1000), 100, 10)
# fm op fm
stopifnot(sum(fm.mapply2(m, m2, fm.bo.add)
			  == fm.mapply2(m, m2, "+")) == length(m))
v <- fm.runif(1000)
v2 <- fm.runif(1000)
# fmV op fmV
stopifnot(sum(fm.mapply2(v, v2, fm.bo.add)
			  == fm.mapply2(v, v2, "+")) == length(v))
vcol <- fm.runif(nrow(m))
# fm op fmV
stopifnot(sum(fm.conv.FM2R(fm.mapply2(m, vcol, fm.bo.add))
			  == fm.conv.FM2R(m) + fm.conv.FM2R(vcol)) == length(m))
mtmp <- matrix(runif(length(m)), nrow(m), ncol(m))
# fm op matrix
stopifnot(sum(fm.conv.FM2R(fm.mapply2(m, mtmp, fm.bo.add))
			  == fm.conv.FM2R(m) + mtmp) == length(m))
# matrix op fm
stopifnot(sum(fm.conv.FM2R(fm.mapply2(mtmp, m, fm.bo.sub))
			  == mtmp - fm.conv.FM2R(m)) == length(m))
vcol <- runif(nrow(m))
# fm op ANY
stopifnot(sum(fm.conv.FM2R(fm.mapply2(m, vcol, fm.bo.add))
			  == fm.conv.FM2R(m) + vcol) == length(m))
stopifnot(sum(fm.conv.FM2R(fm.mapply2(m, 1, fm.bo.add))
			  == fm.conv.FM2R(m) + 1) == length(m))
# ANY op fm
stopifnot(sum(fm.conv.FM2R(fm.mapply2(1, m, fm.bo.sub))
			  == 1 - fm.conv.FM2R(m)) == length(m))
vtmp <- runif(length(v))
# fmV op ANY
stopifnot(sum(fm.conv.FM2R(fm.mapply2(v, vtmp, fm.bo.add))
			  == fm.conv.FM2R(v) + vtmp) == length(v))
stopifnot(sum(fm.conv.FM2R(fm.mapply2(v, 1, fm.bo.add))
			  == fm.conv.FM2R(v) + 1) == length(v))
stopifnot(sum(fm.conv.FM2R(fm.mapply2(vtmp, v, fm.bo.sub))
			  == vtmp - fm.conv.FM2R(v)) == length(v))
stopifnot(sum(fm.conv.FM2R(fm.mapply2(1, v, fm.bo.sub))
			  == 1 - fm.conv.FM2R(v)) == length(v))

# test ifelse
v1 <- fm.runif(1000)
test <- v1 > 0.5
res <- ifelse(test, 0, v1)
rres <- ifelse(fm.conv.FM2R(test), 0, fm.conv.FM2R(v1))
stopifnot(sum(rres == res) == length(v1))

res <- ifelse(test, v1, 0)
rres <- ifelse(fm.conv.FM2R(test), fm.conv.FM2R(v1), 0)
stopifnot(sum(rres == res) == length(v1))

m1 <- fm.matrix(v1, 100, 10)
test <- m1 > 0.5
res <- ifelse(test, 0, m1)
rres <- ifelse(fm.conv.FM2R(test), 0, fm.conv.FM2R(m1))
stopifnot(sum(rres == res) == length(m1))

res <- ifelse(test, m1, 0)
rres <- ifelse(fm.conv.FM2R(test), fm.conv.FM2R(m1), 0)
stopifnot(sum(rres == res) == length(m1))

# Test is.na
v1 <- runif(1000)
v1.int <- as.integer(v1 * 100)
v1[v1 < 0.2] <- NA
v1.int[v1.int < 20] <- NA
cat("#NA in v1:", sum(is.na(v1)), "\n")
cat("#NA in v1.int:", sum(is.na(v1.int)), "\n")
stopifnot(sum(is.na(v1) == is.na(fm.conv.R2FM(v1))) == length(v1))
stopifnot(sum(is.na(v1.int) == is.na(fm.conv.R2FM(v1.int))) == length(v1.int))

m1 <- matrix(v1, 100, 10)
m1.int <- matrix(v1.int, 100, 10)
cat("#NA in m1:", sum(is.na(m1)), "\n")
cat("#NA in m1.int:", sum(is.na(m1.int)), "\n")
stopifnot(sum(is.na(m1) == is.na(fm.conv.R2FM(m1))) == length(m1))
stopifnot(sum(is.na(m1.int) == is.na(fm.conv.R2FM(m1.int))) == length(m1.int))

# Test is.nan
v1 <- runif(1000)
v1.int <- as.integer(v1 * 100)
v1[v1 < 0.2] <- NaN
cat("#NaN in v1:", sum(is.nan(v1)), "\n")
stopifnot(sum(is.nan(v1) == is.nan(fm.conv.R2FM(v1))) == length(v1))
stopifnot(sum(is.nan(v1.int) == is.nan(fm.conv.R2FM(v1.int))) == length(v1.int))

m1 <- matrix(v1, 100, 10)
m1.int <- matrix(v1.int, 100, 10)
cat("#NaN in m1:", sum(is.nan(m1)), "\n")
stopifnot(sum(is.nan(m1) == is.nan(fm.conv.R2FM(m1))) == length(m1))
stopifnot(sum(is.nan(m1.int) == is.nan(fm.conv.R2FM(m1.int))) == length(m1.int))

# Test is.finite and is.infinite
v1 <- runif(1000)
v1.int <- as.integer(v1 * 100)
v1[v1 < 0.2] <- Inf
cat("#Inf in v1:", sum(is.infinite(v1)), "\n")
stopifnot(sum(is.infinite(v1) == is.infinite(fm.conv.R2FM(v1))) == length(v1))
stopifnot(sum(is.infinite(v1.int) == is.infinite(fm.conv.R2FM(v1.int))) == length(v1.int))
stopifnot(sum(is.finite(v1) == is.finite(fm.conv.R2FM(v1))) == length(v1))
stopifnot(sum(is.finite(v1.int) == is.finite(fm.conv.R2FM(v1.int))) == length(v1.int))

m1 <- matrix(v1, 100, 10)
m1.int <- matrix(v1.int, 100, 10)
cat("#Inf in m1:", sum(is.infinite(m1)), "\n")
stopifnot(sum(is.infinite(m1) == is.infinite(fm.conv.R2FM(m1))) == length(m1))
stopifnot(sum(is.infinite(m1.int) == is.infinite(fm.conv.R2FM(m1.int))) == length(m1.int))
stopifnot(sum(is.finite(m1) == is.finite(fm.conv.R2FM(m1))) == length(m1))
stopifnot(sum(is.finite(m1.int) == is.finite(fm.conv.R2FM(m1.int))) == length(m1.int))

# test materialize a matrix
m1 <- fm.matrix(fm.runif(2000), 200, 10)
m2 <- m1 < 0.5
m3 <- fm.materialize(m2)
stopifnot(typeof(m3) == "logical")
m4 <- fm.conv.layout(m2, TRUE)
stopifnot(typeof(m4) == "logical")
m4 <- fm.conv.layout(m2, FALSE)
stopifnot(typeof(m4) == "logical")

# Test on a sparse matrix
print("load sparse matrix")
ig <- read.graph("wiki-Vote1.txt", directed=TRUE)
fg <- fg.load.graph("wiki-Vote1.txt", directed=TRUE)
fm <- fm.get.sparse.matrix(fg)
adj.m <- get.adjacency(ig)
stopifnot(!fm.is.sym(fm))
stopifnot(fm.is.sparse(fm))

print("SpMV")
fm.vec <- fm.rep.int(1, dim(fm)[2])
fm.res <- fm.multiply(fm, fm.vec)
vec <- fm.conv.FM2R(fm.vec)
res <- adj.m %*% vec
stopifnot(sum(fm.conv.FM2R(fm.res) == res) == length(res))

# Test SpMM
print("SpMM")
fm.mat <- fm.matrix(fm.rep.int(1, dim(fm)[2] * 20), dim(fm)[2], 20, TRUE)
fm.res <- fm.multiply(fm, fm.mat)
mat <- fm.conv.FM2R(fm.mat)
res <- adj.m %*% mat
stopifnot(sum(fm.conv.FM2R(fm.res) == res) == length(res))

# Test on a very large dense matrix.
# convert a large dense matrix in R to FlashMatrixR.
