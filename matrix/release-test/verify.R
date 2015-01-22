library(FlashGraphR)
#fg.set.conf("flash-graph/conf/run_test.txt")

# Test a FlashMatrixR vector.
print("convert a FlashMatrixR vector to R vector")
fm.vec <- fm.create.vector(2000, 1)
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
fm.mat <- fm.matrix(fm.create.vector(2000, 1), 20, 100)
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

mat <- matrix(1:2000, 20, 100)

# Test a dense matrix times a vector
fm.vec <- fm.create.vector(fm.ncol(fm.mat), 1)
vec <- fm.conv.FM2R(fm.vec)
res <- mat %*% vec

# TODO
#print("column-wise matrix times a vector")
#fm.mat <- fm.conv.R2FM(mat, byrow = FALSE)
#fm.res <- fm.multiply(fm.mat, fm.vec)
#stopifnot(sum(res[, 0] == fm.conv.FM2R(fm.res)) == fm.length(fm.res))

print("Row-wise matrix times a vector")
fm.mat <- fm.conv.R2FM(mat, byrow = TRUE)
fm.res <- fm.multiply(fm.mat, fm.vec)
stopifnot(sum(res[, 1] == fm.conv.FM2R(fm.res)) == fm.length(fm.res))

print("a wide dense matrix times a tall dense matrix")
left.mat <- matrix(1:2000, 20, 100)
right.mat <- matrix(1:2000, 100, 20)
res <- left.mat %*% right.mat

#print("Column-wise matrix times a column-wise matrix")
#fm.left <- fm.conv.R2FM(left.mat, byrow = FALSE)
#fm.right <- fm.conv.R2FM(right.mat, byrow = FALSE)
#fm.res <- fm.multiply(fm.left, fm.right)
#stopifnot(sum(res == fm.conv.FM2R(fm.res)) == fm.nrow(fm.res) * fm.ncol(fm.res))

#print("Column-wise matrix times a row-wise matrix")
#fm.left <- fm.conv.R2FM(left.mat, byrow = FALSE)
#fm.right <- fm.conv.R2FM(right.mat, byrow = TRUE)
#fm.res <- fm.multiply(fm.left, fm.right)
#stopifnot(sum(res == fm.conv.FM2R(fm.res)) == fm.nrow(fm.res) * fm.ncol(fm.res))

print("Row-wise matrix times a column-wise matrix")
fm.left <- fm.conv.R2FM(left.mat, byrow = TRUE)
fm.right <- fm.conv.R2FM(right.mat, byrow = FALSE)
fm.res <- fm.multiply(fm.left, fm.right)
stopifnot(sum(res == fm.conv.FM2R(fm.res)) == fm.nrow(fm.res) * fm.ncol(fm.res))

#print("Row-wise matrix times a row-wise matrix")
#fm.left <- fm.conv.R2FM(left.mat, byrow = TRUE)
#fm.right <- fm.conv.R2FM(right.mat, byrow = TRUE)
#fm.res <- fm.multiply(fm.left, fm.right)
#stopifnot(sum(res == fm.conv.FM2R(fm.res)) == fm.nrow(fm.res) * fm.ncol(fm.res))

print("a tall dense matrix times a wide dense matrix")
left.mat <- matrix(1:2000, 100, 20)
right.mat <- matrix(1:2000, 20, 10)
res <- left.mat %*% right.mat

print("Column-wise matrix times a column-wise matrix")
fm.left <- fm.conv.R2FM(left.mat, byrow = FALSE)
fm.right <- fm.conv.R2FM(right.mat, byrow = FALSE)
fm.res <- fm.multiply(fm.left, fm.right)
stopifnot(sum(res == fm.conv.FM2R(fm.res)) == fm.nrow(fm.res) * fm.ncol(fm.res))

print("Column-wise matrix times a row-wise matrix")
fm.left <- fm.conv.R2FM(left.mat, byrow = FALSE)
fm.right <- fm.conv.R2FM(right.mat, byrow = TRUE)
fm.res <- fm.multiply(fm.left, fm.right)
stopifnot(sum(res == fm.conv.FM2R(fm.res)) == fm.nrow(fm.res) * fm.ncol(fm.res))

print("Row-wise matrix times a column-wise matrix")
fm.left <- fm.conv.R2FM(left.mat, byrow = TRUE)
fm.right <- fm.conv.R2FM(right.mat, byrow = FALSE)
fm.res <- fm.multiply(fm.left, fm.right)
stopifnot(sum(res == fm.conv.FM2R(fm.res)) == fm.nrow(fm.res) * fm.ncol(fm.res))

#print("Row-wise matrix times a row-wise matrix")
#fm.left <- fm.conv.R2FM(left.mat, byrow = TRUE)
#fm.right <- fm.conv.R2FM(right.mat, byrow = TRUE)
#fm.res <- fm.multiply(fm.left, fm.right)
#stopifnot(sum(res == fm.conv.FM2R(fm.res)) == fm.nrow(fm.res) * fm.ncol(fm.res))

# Test on a sparse matrix
print("load sparse matrix")
fg <- fg.load.graph("data/wiki.adj-v4", "data/wiki.index-v4")
fm <- fm.get.matrix(fg)
stopifnot(!fm.is.sym(fm))
stopifnot(fm.is.sparse(fm))

print("SpMV")
fm.vec <- fm.create.vector(fm.ncol(fm), 1)
fm.res <- fm.multiply(fm, fm.vec)
stopifnot(fm.sum(fm.res) == fg.ecount(fg))

# Test SpMM
# TODO

# Test on float-point vector/matrix

# Test on a very large dense matrix.
# convert a large dense matrix in R to FlashMatrixR.
