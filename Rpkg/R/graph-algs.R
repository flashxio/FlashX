PageRank <- function(A, d, max.niters) {
	N <- dim(A)[1]
	epsilon <- 1e-2/N
	cat("There are", N, "vertices. eps =", epsilon)
	one <- fm.rep.int(1, N)
	out.deg <- A %*% one
	one <- NULL
	pr1 <- fm.rep.int(1/N, N)
	converge <- 0

	niters <- 0
	A <- t(A)
	while (converge < N && niters < max.niters) {
		pr2 <- (1-d)/N+d*(A %*% (pr1/out.deg))
		diff <- abs(pr1-pr2)
		converge <- sum(diff < epsilon)
		pr1 <- pr2
		niters <- niters + 1
	}
	pr2
}
