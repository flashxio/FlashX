PageRank <- function(A, d = 0.15, max.niters = 30, epsilon = 1e-2,
					 verbose = FALSE)
{
	orig.test.na <- fm.env$fm.test.na
	fm.set.test.na(FALSE)
	N <- dim(A)[1]
	epsilon <- epsilon / N
	cat("There are", N, "vertices.", "\n")
	one <- fm.rep.int(1, N)
	out.deg <- A %*% one
	one <- NULL
	pr1 <- fm.rep.int(1/N, N)
	L1 <- N

	niters <- 0
	A <- t(A)
	while (L1 >= epsilon && niters < max.niters) {
		pr2 <- d/N+(1-d)*(A %*% (pr1/out.deg))
		gc()
		pr2 <- fm.materialize(pr2)
		diff <- abs(pr1-pr2)
		L1 <- sum(diff)
		if (verbose)
			cat("iter", niters, ", L1:", L1 * N, "\n")
		pr1 <- pr2
		niters <- niters + 1
	}
	fm.set.test.na(orig.test.na)
	pr2
}
