# This file contains some implementations from the mvtnorm package.

fm.dmvnorm <- function(X, mu, covar, log=FALSE)
{
	orig.test.na <- .env.int$fm.test.na
	.set.test.na(FALSE)
	if (fm.is.matrix(covar))
		covar <- fm.conv.FM2R(covar)
	covar.inv <- solve(covar)
	X1 <- sweep(X, 2, mu, "-")
	X2 <- X1 %*% covar.inv
	X3 <- fm.agg.mat(X2 * X1, 1, "+")
	k <- dim(covar)[1]
	dec <- tryCatch(chol(covar), error = function(e) e)
	logdet <- -sum(log(diag(dec)))
	logret <- logdet - 0.5 * k * log(2 * pi) - 0.5 * X3
	if (log)
		ret <- logret
	else
		ret <- exp(logret)
	.set.test.na(orig.test.na)
	ret
}
