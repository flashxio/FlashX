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

# This file contains some implementations from the mvtnorm package.

fm.dmvnorm <- function(X, mu, covar, log=FALSE)
{
	if (fm.is.matrix(covar))
		covar <- fm.conv.FM2R(covar)
	covar.inv <- solve(covar)
	X1 <- sweep(X, 2, mu, "-")
	X2 <- X1 %*% covar.inv
	X3 <- fm.agg.mat(X2 * X1, 1, "+")
	k <- dim(covar)[1]
	dec <- tryCatch(chol(covar), error = function(e) e)
	logdet <- -sum(log(diag(dec)))
	logret <- logdet - 0.5 * k * log(2 * pi) - 0.5 * X3
	if (log)
		ret <- logret
	else
		ret <- exp(logret)
	ret
}
