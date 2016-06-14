logistic.grad <- function(x, y, theta)
{
	gradient <- (1 / length(y)) * (t(x) %*% (1/(1 + exp(-x %*% t(theta))) - y))
	return(t(fm.conv.FM2R(gradient)))
}

logistic.cost <- function(x, y, theta)
{
	m <- length(y)
	(1/m)*sum(y*(-x %*% t(theta)) + log(1 + exp(x %*% t(theta))))
}

logistic.regression <- function(x, y)
{
	gradient.descent(x, y, logistic.grad, cost=logistic.cost, alpha=0.001)
}
