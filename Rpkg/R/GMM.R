comp.prob <- function(X, mus, covars, phi)
{
	k <- length(covars)

	log.likely.list <- list()
	for (i in 1:k)
		log.likely.list[[i]] <- fm.dmvnorm(X, mus[,i], covars[[i]], TRUE)
	log.likely <- fm.cbind.list(log.likely.list)
	# I need to materialize a matrix here first to speed up. Why?
	log.likely <- fm.materialize(log.likely)

	P.list <- list()
	for (i in 1:k)
		P.list[[i]] <- phi[i] / (exp(log.likely - log.likely[,i]) %*% phi)
	P <- fm.cbind.list(P.list)
	P <- fm.materialize(P)

	max.log.likely <- fm.agg.mat(log.likely, 1, fm.bo.max)
	rel.log.likely <- log.likely - max.log.likely
	rel.likely <- exp(rel.log.likely)
	log.like <- sum(log(rel.likely %*% phi)) + sum(max.log.likely)

	list(P=P, log.like=log.like)
}

GMM <- function(X, k, maxiters, verbose=FALSE)
{
	orig.test.na <- fm.env$fm.test.na
	fm.set.test.na(FALSE)

	m <- dim(X)[1]

	# Random init
	# TODO alternatively, we can use KMeans to initialize it.
	rand.k <- as.integer(runif(k, 1, m))
	mus <- t(X[rand.k,])
	init.covar <- cov(X)
	covars <- list()
	for (i in 1:k)
		covars[[i]] <- fm.conv.FM2R(init.covar)
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
		mus <- sweep(t(X) %*% P, 2, phi * m, "/")
		for (j in 1:k)
			covars[[j]] <- fm.conv.FM2R(cov.wt(X, P[,j])$cov)
		end.t <- Sys.time()
		if (verbose)
			cat("M-step takes", as.integer(end.t) - as.integer(start.t),
				"seconds\n")
	}
	mus <- fm.conv.FM2R(mus)
	clust.ids <- fm.materialize(fm.agg.mat(P, 1, fm.bo.which.max))
	fm.set.test.na(orig.test.na)
	list(P=P, clust.ids=clust.ids, niters=iter, phi=phi, mus=mus, covars=covars)
}
