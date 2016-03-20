GMM <- function(X, k, maxiters, verbose=FALSE)
{
	log.like <- function(X, phi, mus, covars) {
		P.list <- list()
		phi <- fm.conv.FM2R(phi)
		for (i in 1:length(covars))
			P.list[[i]] <- phi[i] * fm.dmvnorm(X, mus[,i], covars[[i]])
		P <- fm.cbind.list(P.list)
		sum(log(rowSums(P)))
	}

	m <- dim(X)[1]

	# Random init
	# TODO alternatively, we can use KMeans to initialize it.
	rand.k <- as.integer(runif(k, 1, m))
	mus <- t(X[rand.k,])
	init.covar <- cov(X)
	covars <- list()
	for (i in 1:k)
		covars[[i]] <- init.covar
	phi <- fm.rep.int(1/m, k)

	old.state <- list()
	for (iter in 1:maxiters) {
		# E-step
		P.list <- list()
		for (i in 1:k) {
			curr.covar <- covars[[i]]
			P.list[[i]] <- fm.dmvnorm(X, mus[,i], covars[[i]])
		}
		P <- fm.cbind.list(P.list)
		P <- sweep(P, 2, phi, "*") / (P %*% phi)

		# M-step
		phi <- colSums(P)/m
		mus <- sweep(t(X) %*% P, 2, phi * m, "/")
		for (j in 1:k)
			covars[[j]] <- cov.wt(X, P[,j])$cov

		if (length(old.state) > 0) {
			new.like <- log.like(X, phi, mus, covars)
			old.like <- log.like(X, old.state$phi, old.state$mus, old.state$covars)
			if (verbose)
				cat("log likelihood:", new.like, "\n")
			if (new.like - old.like < 0.01)
				break
		}
		old.state <- list(phi=phi, mus=mus, covars=covars)
	}
	P
}
