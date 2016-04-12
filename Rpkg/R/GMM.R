gmm.covs <- function (x, wts)
{
	k <- ncol(wts)
	n <- nrow(x)
	if (nrow(wts) != n)
		stop("length of 'wt' must equal the number of rows in 'x'")

	s <- fm.colSums(wts, TRUE)
	s2 <- fm.colSums(wts * wts, TRUE)
	ret <- fm.materialize(s, s2)
	s <- fm.conv.FM2R(ret[[1]])
	s2 <- fm.conv.FM2R(ret[[2]])

	centers <- list()
	xs.cp <- list()
	for (i in 1:k) {
		wx <- wts[,i] * x / s[i]
		# fm.colSums(wt * x, TRUE)
		centers[[i]] <- fm.rowSums(t(wx), TRUE)
		xs.cp[[i]] <- fm.crossprod(wx, x, lazy=TRUE)
	}
	centers <- fm.materialize.list(centers)
	xs.cp <- fm.materialize.list(xs.cp)

	for (i in 1:k) {
		center <- fm.conv.FM2R(centers[[i]])
		x.cp <- fm.conv.FM2R(xs.cp[[i]])
		x.cp <- x.cp - center %*% t(center)
		# x.cp/(1 - sum(wt^2))
		xs.cp[[i]] <- x.cp/(1 - s2[i]/s[i]/s[i])
	}
	xs.cp
}

comp.prob <- function(X, mus, covars, phi)
{
	k <- length(covars)

	log.likely.list <- list()
	for (i in 1:k)
		log.likely.list[[i]] <- fm.dmvnorm(X, mus[,i], covars[[i]], TRUE)
	log.likely <- fm.cbind.list(log.likely.list)
	# I need to materialize a matrix here first to speed up. Why?
	log.likely <- fm.materialize(log.likely)

	max.log.likely <- fm.agg.mat(log.likely, 1, fm.bo.max)
	rel.likely <- exp(log.likely - max.log.likely) %*% phi
	sum1 <- fm.sum(log(rel.likely), TRUE)
	sum2 <- fm.sum(max.log.likely, TRUE)

	P <- sweep(exp(log.likely - max.log.likely), 2, phi, "*") / rel.likely
	ret <- fm.materialize(P, sum1, sum2)
	log.like <- fm.conv.FM2R(ret[[2]]) + fm.conv.FM2R(ret[[3]])

	list(P=ret[[1]], log.like=log.like)
}

GMM <- function(X, k, maxiters, verbose=FALSE)
{
	orig.test.na <- fm.env$fm.test.na
	fm.set.test.na(FALSE)

	m <- dim(X)[1]

	X <- fm.conv.layout(X, FALSE)
	X <- fm.materialize(X)

	if (!all(is.finite(X)))
		stop("'x' must contain finite values only")

	# Random init
	# TODO alternatively, we can use KMeans to initialize it.
	rand.k <- as.integer(runif(k, 1, m))
	mus <- fm.conv.FM2R(t(X[rand.k,]))
	init.covar <- cov(X)
	covars <- list()
	for (i in 1:k)
		covars[[i]] <- init.covar
	phi <- rep.int(1/m, k)

	for (iter in 1:maxiters) {
		if (verbose)
			cat("iter", iter, "\n")
		# E-step
		start.t <- Sys.time()
		ret <- comp.prob(X, mus, covars, phi)
		P <- ret$P
		new.like <- ret$log.like
		if (iter > 1) {
			if (verbose) {
				cat("new log likelihood:", new.like, ", old log likelihood:",
					old.like, "\n")
			}
			if (new.like - old.like < 0.01)
				break
		}
		old.like <- new.like
		gc()
		end.t <- Sys.time()
		if (verbose)
			cat("E-step takes", as.integer(end.t) - as.integer(start.t),
				"seconds\n")

		# M-step
		start.t <- Sys.time()
		phi <- fm.conv.FM2R(colSums(P)/m)
		if (verbose)
			print(phi)
		mus <- sweep(fm.conv.FM2R(t(X) %*% P), 2, phi * m, "/")

		covars <- gmm.covs(X, P)
		end.t <- Sys.time()
		if (verbose)
			cat("M-step takes", as.integer(end.t) - as.integer(start.t),
				"seconds\n")
	}
	clust.ids <- fm.materialize(fm.agg.mat(P, 1, fm.bo.which.max))
	fm.set.test.na(orig.test.na)
	list(P=P, clust.ids=clust.ids, niters=iter, phi=phi, mus=mus, covars=covars)
}
