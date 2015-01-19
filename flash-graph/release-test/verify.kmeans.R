require(FlashGraphR)
require(MASS)
require(mclust)

plot.color.func <- function() {
	# TODO
}

test.kmeans <- function(k, n=50, std.dev=.3, min.ari.thresh=.5)
{
	# Make data
	mu.vec <- mvrnorm(k, c(0,0), std.dev*diag(2)); # Cluster means
	sigma <- diag(2)

	data <- mvrnorm(n, mu.vec[1,], sigma)
	print(paste("Set a cluster mean to: ", toString(mu.vec[1,])))
	truth <- c(replicate(n, 1)) # True cluster assignments

	for (i in 2:k) {
		data <- rbind(data, mvrnorm(n, mu.vec[i,], sigma))
		truth <- c(truth, replicate(n, i)) # These are the true clusters
	}

	fg.kms <- fg.kmeans(data, k, max.iters=100)
	R.kms <- kmeans(data, k, iter.max=100, algorithm="Forgy")

	#browser()
	fg.ari <- adjustedRandIndex(fg.kms$cluster, truth)
	R.ari <- adjustedRandIndex(R.kms$cluster, truth)

	# Organize the data
	print(paste("The ARI of fg is", fg.ari))
	print(paste("The ARI of R kms is", R.ari))

	# Plot result

	#plot(x, data[,1], data[,2], )

	stopifnot(fg.ari >= min.ari.thresh)
}


test.kmeans(5)
