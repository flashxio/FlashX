# Copyright 2016 Open Connectome Project (http://openconnecto.me)
# Written by Da Zheng (zhengda1936@gmail.com)
#
# This file is part of FlashR.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

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

logistic.regression <- function(X, y, method=c("GD", "Newton", "LS", "RNS", "Uniform"),
								hessian_size=0.1)
{
	if (method == "GD")
		get.hessian <- NULL
	else
		get.hessian <- logistic.hessian
	params <- list(c=0.5, ro=0.2, linesearch=is.null(get.hessian),
				   num.iters=500, out.path=FALSE, method=method, hessian_size=hessian_size)
	gradient.descent(X, y, logistic.grad, get.hessian, cost=logistic.cost, params)
}
