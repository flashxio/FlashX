gradient.descent <- function(x, y, grad, cost=NULL, alpha=0.1,
							 num.iterations=500, threshold=1e-5,
							 output.path=FALSE)
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
		g <- grad(x, y, theta)
		cat("L2 of gradient:", L2(g), "\n")
		theta <- theta - alpha * g
		if(all(is.na(theta))) break
		theta.path <- rbind(theta.path, theta)
		if (!is.null(cost))
			cat("cost:", cost(x, y, theta), "\n")
		if(i > 2) if(all(abs(theta - theta.path[i-1,]) < threshold)) break 
	}

	if(output.path) return(theta.path) else return(theta.path[nrow(theta.path),])
}
