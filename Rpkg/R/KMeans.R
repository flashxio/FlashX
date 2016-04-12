fm.KMeans <- function(data, K, max.iters=10, debug=FALSE)
{
	orig.test.na <- fm.env$fm.test.na
	fm.set.test.na(FALSE)

	n <- dim(data)[1]
	m <- dim(data)[2]
	agg.sum <- fm.create.agg.op(fm.bo.add, fm.bo.add, "sum")
	agg.which.min <- fm.create.agg.op(fm.bo.which.min, NULL, "which.min")

	cal.centers <- function(data, parts) {
		centers1 <- fm.groupby(data, 2, parts, agg.sum)
		centers1 <- as.matrix(centers1)
		cnts <- fm.table(parts)
		centers.idx <- as.vector(cnts$val) + 1
		# Some centers may not have been initialized if no data points
		# are assigned to them. I need to reset those data centers.
		centers <- matrix(rep.int(0, length(centers1)), dim(centers1)[1],
						  dim(centers1)[2])
		centers[centers.idx,] <- centers1[centers.idx,]
		# Calculate the mean for each cluster, which is the cluster center
		# for the next iteration.
		sizes <- rep.int(1, dim(centers)[1])
		sizes[centers.idx] <- as.vector(cnts$Freq)
		centers <- diag(1/sizes) %*% centers
		fm.as.matrix(centers)
	}

#	parts <- fm.as.integer(floor(fm.runif(n, min=0, max=K)))
#	new.centers <- cal.centers(data, fm.as.factor(parts, K))
	rand.k <- runif(K, 1, nrow(data))
	new.centers <- data[rand.k,]
	parts <- NULL

	iter <- 0
	start.time <- Sys.time()
	num.moves <- nrow(data)
	while (num.moves > 0 && iter < max.iters) {
		if (debug)
			iter.start <- Sys.time()
		centers <- new.centers
		old.parts <- parts
		gc()
		m <- fm.inner.prod(data, t(centers), fm.bo.euclidean, fm.bo.add)
		parts <- fm.as.integer(fm.agg.mat(m, 1, agg.which.min) - 1)
		# Have the vector materialized during the computation.
		fm.set.materialize.level(parts, 2, TRUE)

		new.centers <- cal.centers(data, fm.as.factor(parts, K))
		if (!is.null(old.parts))
			num.moves <- sum(fm.as.numeric(old.parts != parts))
		iter <- iter + 1
		if (debug) {
			iter.end <- Sys.time()
			cat("iteration", iter, "takes",
				as.numeric(iter.end) - as.numeric(iter.start),
				"seconds and moves", num.moves, "data points\n")
		}
		old.parts <- NULL
		gc()
	}
	end.time <- Sys.time()
	cat("KMeans takes", iter , "iterations and",
		as.numeric(end.time) - as.numeric(start.time), "seconds\n")
	fm.set.test.na(orig.test.na)
	parts
}
