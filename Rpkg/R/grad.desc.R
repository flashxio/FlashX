gradient.descent <- function(X, y, get.grad, get.hessian, cost, params)
{
	# Add x_0 = 1 as the first column
	m <- if(is.vector(X)) length(X) else nrow(X)
	if(is.vector(X) || (!all(X[,1] == 1))) X <- fm.cbind(fm.rep.int(1, m), X)

	num.features <- ncol(X)

	# Initialize the parameters
	theta <- matrix(rep(0, num.features), nrow=1)

	L2 <- function(X) {
		sqrt(sum(X * X))
	}

	# Look at the values over each iteration
	theta.path <- theta
	for (i in 1:params$num.iters) {
		g <- get.grad(X, y, theta)
		if (!is.null(get.hessian) && i > 5) {
			D2 <- get.hessian(X, y, theta)
			H <- fm.conv.FM2R(t(X) %*% sweep(X, 1, fm.as.vector(D2), FUN="*"))
			H <- H / min(abs(H))
			z <- pcg(H, as.vector(-g))
			z <- t(z)
			params$linesearch <- FALSE
		}
		else {
			params$linesearch <- TRUE
			z <- -g
		}

		eta <- 1
		if (params$linesearch) {
			eta <- 0.01
			l <- cost(X, y, theta)
			delta = params$c * t(z) %*% g
			while (cost(X, y, theta + eta * z) >= l + delta * eta)
				eta <- eta * params$ro
			print(eta)
		}
		if (eta == 1)
			params$linesearch <- FALSE

		theta <- theta + z * eta
		if(all(is.na(theta))) break
		theta.path <- rbind(theta.path, theta)
		cat(i,  ": L2(g) =", L2(g), ", cost:", cost(X, y, theta), "\n")
	}

	if(params$out.path) return(theta.path) else return(theta.path[nrow(theta.path),])
}
