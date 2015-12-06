FG <- TRUE
verbose <- FALSE
fgcheck <- function() {
	if ("FlashR" %in% rownames(installed.packages()) == FALSE) {
		assign("FG", FALSE, envir = .GlobalEnv)
	} else {
		require(FlashR)
	}
}

fgcheck()
require(MASS)
require(mclust)
ceil <- ceiling # Funcdef

# Plot with a colorscheme that is dictated by which cluster each
# point falls into
plot.by.cluster <- function (data, clusters, title, xlims, ylims, colors, k) {
	plot.called <- FALSE
	dev.new()

	for (i in 1:k) {
		if (!plot.called) {
			plot(data[which(clusters == i), 1], data[which(clusters == i),2], 
				 main=title,col=colors[i], xlim=xlims, ylim=ylims)
			plot.called <- TRUE
		} else {
			points(data[which(clusters == i), 1], data[which(clusters == i),2],
				   col=colors[i], xlim=xlims, ylim=ylims)
		}
	}
}

# Allow no more than 10% degredation compared to R.kmeans
test.kmeans <- function(k, n=50, std.dev=10, min.ari.thresh=.1)
{
	### Make data ####
	mu.vec <- mvrnorm(k, c(0,0), std.dev*diag(2)); # Cluster means

	if (verbose) {
		print(paste("The sum of variance btwn mean dims is", 
					toString(var(mu.vec[,1]) + var(mu.vec[,2]))))
		print(paste("Means variance:", toString(var(mu.vec))))
	}

	sigma <- diag(2)

	colors <- rainbow(k)
	data <- mvrnorm(n, mu.vec[1,], sigma) # class 1 distibution

	if (verbose) {
		print(paste("Set a cluster mean to: ", toString(mu.vec[1,])))
		cat("\n")
	}

	truth <- c(replicate(n, 1)) # True cluster assignments

	if (verbose) {
		dev.new()
		# Plot the first points from class 1
		plot(data[,1], data[,2], col =colors[1], main="Truth", 
			 xlim=c(-std.dev, std.dev), ylim=c(-std.dev,std.dev))

		print(paste("Sum of variance dims, cluster 1: ", toString(var(data[,1]) + var(data[,2]))))
		print(paste("Cluster 1 var:", toString(var(data))))
		cat("\n")
	}

	for (i in 2:k) {
		data <- rbind(data, mvrnorm(n, mu.vec[i,], sigma))
		truth <- c(truth, replicate(n, i)) # These are the true clusters

		if (verbose) {
			points(data[(((i-1)*n)+1):dim(data)[1],1], data[(((i-1)*n)+1):dim(data)[1],2], 
				   col=colors[i]) 

			print(paste("Sum of variance dims, cluster", toString(i),
						":",toString( var(data[(((i-1)*n)+1):dim(data)[1],1]) + 
									 var(data[(((i-1)*n)+1):dim(data)[1],2]) )))
			print(paste("Cluster", toString(i), "var:", 
						toString(var(data[(((i-1)*n)+1):dim(data)[1],]))))
			cat("\n")
			print(paste("Set a cluster mean to: ", toString(mu.vec[i,])))
		}
	}

	R.kms <- kmeans(data, k, iter.max=100, algorithm="Forgy")
	R.ari <- adjustedRandIndex(R.kms$cluster, truth)
	print(paste("The ARI of R kms is", R.ari))

	if (FG) {
		fg.kms <- fg.kmeans(data, k, max.iters=100)
		fg.ari <- adjustedRandIndex(fg.kms$cluster, truth)
		print(paste("The ARI of fg is", fg.ari))
		stopifnot((fg.ari > R.ari) || 
				  (abs(fg.ari - R.ari)/mean(c(fg.ari,R.ari))) < min.ari.thresh) # 10% difference only
	}

	if (verbose) {
		# Cause we have all the data we get the limits
		xlims <- c(floor(min(data[,1])), ceil(max(data[,1]))) 
		ylims <- c(floor(min(data[,2])), ceil(max(data[,2])))

		# Plot result
		plot.by.cluster(data, R.kms$cluster, "R-clusters colored by cluster",
						xlims, ylims, colors, k)
		if (FG) {
			plot.by.cluster(data, fg.kms$cluster, "FG-clusters colored by cluster", 
							xlims, ylims, colors, k)
		}
	}
}

# Iris is a dataset within R
test.with.iris <- function() {
	data <- data.matrix(iris[,1:dim(iris)[2]-1])
	colnames(data) <- NULL

	# labels
	# setosa = 1
	# versicolor = 2
	# virginica = 3

	truth <- c(rep(1,50), rep(2,50), rep(3,50))

	R.kms <- kmeans(data, 3, iter.max=100)
	fg.kms <- fg.kmeans(data, 3, max.iters=100)

	print(paste("The ARI of R.kms = ", adjustedRandIndex(R.kms$cluster, truth)))
	print(paste("The ARI of fg.kms =", adjustedRandIndex(fg.kms$cluster, truth)))
}
