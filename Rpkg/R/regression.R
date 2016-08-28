logistic.grad <- function(x, y, theta)
{
	gradient <- (1 / length(y)) * (t(x) %*% (1/(1 + exp(-x %*% t(theta))) - y))
	return(t(fm.conv.FM2R(gradient)))
}

logistic.hessian <- function(x, y, w)
{
	exw <- exp(x %*% t(w))
	ifelse(is.infinite(exw), 0, exw/((1+exw)^2))
}

logistic.cost <- function(x, y, theta)
{
	m <- length(y)
	(1/m)*sum(y*(-x %*% t(theta)) + log(1 + exp(x %*% t(theta))))
}

logistic.regression <- function(x, y, method=c("GD", "Newton"))
{
	if (method == "GD")
		get.hessian <- NULL
	else
		get.hessian <- logistic.hessian
	gradient.descent(x, y, logistic.grad, get.hessian,
					 cost=logistic.cost, alpha=0.00000001)
}
