gradient.descent <- function(x, y, get.grad, get.hessian, cost=NULL, alpha=0.1,
							 num.iterations=500, output.path=FALSE)
{
	# Add x_0 = 1 as the first column
	m <- if(is.vector(x)) length(x) else nrow(x)
	if(is.vector(x) || (!all(x[,1] == 1))) x <- fm.cbind(fm.rep.int(1, m), x)

	num.features <- ncol(x)

	# Initialize the parameters
	theta <- matrix(rep(0, num.features), nrow=1)

	L2 <- function(x) {
		sqrt(sum(x * x))
	}

	# Look at the values over each iteration
	theta.path <- theta
	for (i in 1:num.iterations) {
		g <- get.grad(x, y, theta)
		if (!is.null(get.hessian) && i > 5) {
			D2 <- get.hessian(x, y, theta)
			H <- fm.conv.FM2R(t(x) %*% sweep(x, 1, fm.as.vector(D2), FUN="*"))
			H <- H / min(abs(H))
			z <- pcg(H, as.vector(-g))
			z <- t(z)
		}
		else
			z <- -alpha * g

		# TODO we might need line search here.
		theta <- theta + z
		if(all(is.na(theta))) break
		theta.path <- rbind(theta.path, theta)
		if (is.null(cost))
			cat(i,  ": L2(g) =", L2(g), "\n")
		else
			cat(i,  ": L2(g) =", L2(g), ", cost:", cost(x, y, theta), "\n")
	}

	if(output.path) return(theta.path) else return(theta.path[nrow(theta.path),])
}
