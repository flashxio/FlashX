cal.centers <- function(data, parts) {
	centers1 <- fm.groupby(data, 2, parts, agg.sum)
	centers1 <- as.matrix(centers1)
	cnts <- fm.table(parts)
	centers.idx <- as.vector(cnts$val) + 1
	# Some centers may not have been initialized if no data points are assigned to them.
	# I need to reset those data centers.
	centers <- matrix(rep.int(0, length(centers1)), dim(centers1)[1], dim(centers1)[2])
	centers[centers.idx,] <- centers1[centers.idx,]
	# Calculate the mean for each cluster, which is the cluster center for the next iteration.
	sizes <- rep.int(1, dim(centers)[1])
	sizes[centers.idx] <- as.vector(cnts$Freq)
	centers <- diag(1/sizes) %*% centers
	fm.as.matrix(centers)
}

KMeans <- function(data, K)
{
	data <- fm.conv.layout(data, TRUE)
	data <- fm.materialize(data)
	n <- dim(data)[1]
	m <- dim(data)[2]
	agg.sum <- fm.create.agg.op(fm.bo.add, fm.bo.add, "sum")
	agg.which.min <- fm.create.agg.op(fm.bo.which.min, NULL, "which.min")
	parts <- fm.as.integer(floor(fm.runif(n, min=0, max=K)))
	new.centers <- cal.centers(data, fm.as.factor(parts, K))
	centers <- fm.matrix(fm.rep.int(0, K * m), K, m)
	old.parts <- fm.rep.int(0, n)

	iter <- 0
	start.time <- Sys.time()
	while (sum(old.parts != parts) > 0) {
		centers <- new.centers
		old.parts <- parts
		m <- fm.inner.prod(data, t(centers), fm.bo.euclidean, fm.bo.add)
		parts <- fm.as.integer(fm.agg.mat(m, 1, agg.which.min) - 1)
		new.centers <- cal.centers(data, fm.as.factor(parts, K))
		iter <- iter + 1
	}
	end.time <- Sys.time()
	cat("KMeans takes", iter , "iterations and", end.time - start.time, "seconds\n")
}
