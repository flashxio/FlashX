logistic.grad <- function(X, y, w)
{
	gradient <- (1 / length(y)) * (t(X) %*% (1/(1 + exp(-X %*% t(w))) - y))
	return(t(fm.conv.FM2R(gradient)))
}

logistic.hessian <- function(X, y, w)
{
	exw <- exp(X %*% t(w))
	ifelse(is.infinite(exw), 0, exw/((1+exw)^2))
}

logistic.cost <- function(X, y, w)
{
	m <- length(y)
	xw <- X %*% t(w)
	(1/m)*sum(y*(-xw) + log(1 + exp(xw)))
}

logistic.regression <- function(X, y, method=c("GD", "Newton"))
{
	if (method == "GD")
		get.hessian <- NULL
	else
		get.hessian <- logistic.hessian
	gradient.descent(X, y, logistic.grad, get.hessian,
					 cost=logistic.cost, alpha=0.00000001)
}
