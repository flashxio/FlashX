require(igraph)
library("FlashGraphR")
library("irlba")

ig.local.scan <- function(g, order)
{
	interval = 50
	ret <- rep(0, vcount(g))
	start <- 1
	while (start <= vcount(g)) {
		end <- min(start + interval - 1, vcount(g))
		ret[start:end] <- sapply(graph.neighborhood(g, order, nodes=start:end,
													mode="all"), ecount)
		start <- start + interval
	}
	ret
}

verify.cc <- function(fg.res, ig.res)
{
	# get unique components IDs
	uniq.fg.res <- unique(fg.res)
	uniq.ig.res <- unique(ig.res)

	# get the number of isolates from the FlashGraph results.
	fg.no.iso <- sum(fg.res == -1)

	# get the number of unique component IDs
	ig.uniq.no <- length(uniq.ig.res)
	fg.uniq.no <- 0
	if (fg.no.iso == 0) {
		fg.uniq.no <- length(uniq.fg.res)
	} else {
		fg.uniq.no <- length(uniq.fg.res) - 1 + fg.no.iso
	}
	stopifnot(fg.uniq.no == ig.uniq.no)

	# counts the size of each components from FlashGraph
	fg.res.counts <- table(fg.res)
	if (fg.no.iso == 0) {
		fg.uniq.no <- length(fg.res.counts)
		fg.sorted <- sort(fg.res.counts[1:fg.uniq.no])
	}
	else {
		# exclude the isolates from FlashGraph
		fg.uniq.no <- length(fg.res.counts) - 1
		# get the size of each components from FlashGraph, excluding isolates
		fg.sorted <- sort(fg.res.counts[(1:fg.uniq.no) + 1])
	}
	# get the size of each components from iGraph, excluding isolates
	ig.sorted <- sort(table(ig.res))[(ig.uniq.no - fg.uniq.no + 1):ig.uniq.no]
	cmp.res <- fg.sorted == ig.sorted
	stopifnot(sum(cmp.res) == length(cmp.res))
}

check.vectors <- function(name, fg.res, ig.res)
{
	cmp.res <- fg.res == ig.res
	if (sum(cmp.res) != length(cmp.res)) {
		cat(name, "doesn't pass the test\n")
		stopifnot(FALSE)
	}
}

check.matrices <- function(name, fg.res, ig.res)
{
	check.vectors(name, fg.res, ig.res)
}

test.directed <- function(fg, ig)
{
	# test degree
	# this can be used to test the correctness of the generated graph.
	print("test directed degree")
	fg.res <- fg.degree(fg)
	ig.res <- degree(ig)
	check.vectors("degree_test", fg.res, ig.res)
	fg.res <- fg.degree(fg, mode="out")
	ig.res <- degree(ig, mode="out")
	check.vectors("degree_test", fg.res, ig.res)
	fg.res <- fg.degree(fg, mode="in")
	ig.res <- degree(ig, mode="in")
	check.vectors("degree_test", fg.res, ig.res)

	# test coreness
	print("test coreness")
	fg.res <- fg.kcore(fg, 1, 0)
	ig.res <- graph.coreness(ig, mode="all")
	check.vectors("coreness_test", fg.res, ig.res)

	# test WCC
	print("test WCC")
	fg.res <- fg.clusters(fg, mode="weak")
	ig.res <- clusters(ig, mode="weak")$membership
	verify.cc(fg.res, ig.res)

	# test SCC
	print("test SCC")
	fg.res <- fg.clusters(fg, mode="strong")
	ig.res <- clusters(ig, mode="strong")$membership
	verify.cc(fg.res, ig.res)

	# test PageRank
	print("test PageRank")
	fg.res <- fg.page.rank(fg)
	ig.res <- page.rank.old(ig, eps=0.01, old=TRUE)
	num <- sum((abs(fg.res - ig.res) / abs(fg.res)) < 0.02)
	cat("# vertices whose PR diff <= 2% is", num, ", # vertices:", vcount(ig), "\n")

	# test locality scan
	print("test locality statistics")
	fg.res <- fg.local.scan(fg, 2)
	ig.res <- ig.local.scan(ig, 2)
	check.vectors("local-scan2_test", fg.res, ig.res)
	fg.res <- fg.local.scan(fg)
	ig.res <- ig.local.scan(ig, 1)
	check.vectors("local-scan_test", fg.res, ig.res)

	# test topK scan
	print("test topK locality statistics")
	fg.res <- fg.topK.scan(fg, K=10)
	ig.res <- sort(ig.res, decreasing=TRUE)[1:10]
	check.vectors("topK-scan_test", fg.res$scan, ig.res)

	# test matrix multiplication.
	print("A * x");
	x <- runif(vcount(ig), 0, 1)
	ig.matrix <- get.adjacency(ig)
	fg.res <- fg.multiply(fg, x)
	ig.res <- ig.matrix %*% x
	check.vectors("A * x", fg.res, ig.res)

	print("t(A) * x");
	fg.res <- fg.multiply(fg, x, TRUE)
	ig.res <- t(as.matrix(ig.matrix)) %*% x
	check.vectors("t(A) * x", fg.res, ig.res)

	print("A * m");
	x <- matrix(runif(vcount(ig) * 5, 0, 1), nrow=vcount(ig), ncol=5)
	fg.res <- fg.multiply.matrix(fg, x)
	ig.res <- ig.matrix %*% x
	check.matrices("A * m", fg.res, ig.res)

	print("t(A) * m")
	fg.res <- fg.multiply.matrix(fg, x, TRUE)
	ig.res <- t(as.matrix(ig.matrix)) %*% x
	check.matrices("t(A) * m", fg.res, ig.res)
}

test.undirected <- function(fg, ig)
{
	# test degree
	print("test undirected degree")
	fg.res <- fg.degree(fg)
	ig.res <- degree(ig)
	check.vectors("degree_test", fg.res, ig.res)

	# test triangles
	print("test triangle counting on an undirected graph")
	fg.res <- fg.triangles(fg)
	ig.res <- adjacent.triangles(ig)
	check.vectors("undirected-triangle_test", fg.res, ig.res)

	# test locality scan
	print("test locality statistics")
	fg.res <- fg.local.scan(fg)
	ig.res <- sapply(graph.neighborhood(ig, 1, mode="all"), ecount)
	check.vectors("local-scan_test", fg.res, ig.res)

	# test transitivity
	print("test local transitivity")
	fg.res <- fg.transitivity(fg, type="local")
	ig.res <- transitivity(ig, type="local")
	fg.res[is.nan(fg.res)] <- 0
	ig.res[is.nan(ig.res)] <- 0
	check.vectors("local transitivity", fg.res, ig.res)

	print("test global transitivity")
	fg.res <- fg.transitivity(fg, type="global")
	ig.res <- transitivity(ig, type="global")
	stopifnot(fg.res == ig.res)

	# test matrix multiplication.
	print("A * x");
	x <- runif(vcount(ig), 0, 1)
	ig.matrix <- get.adjacency(ig)
	fg.res <- fg.multiply(fg, x)
	ig.res <- ig.matrix %*% x
	check.vectors("A * x", fg.res, ig.res)

	print("t(A) * x");
	fg.res <- fg.multiply(fg, x, TRUE)
	check.vectors("t(A) * x", fg.res, ig.res)

	print("A * m");
	x <- matrix(runif(vcount(ig) * 5, 0, 1), nrow=vcount(ig), ncol=5)
	fg.res <- fg.multiply.matrix(fg, x)
	ig.res <- ig.matrix %*% x
	check.matrices("A * m", fg.res, ig.res)

	print("t(A) * m")
	fg.res <- fg.multiply.matrix(fg, x, TRUE)
	check.matrices("t(A) * m", fg.res, ig.res)
}

test.eigen.directed <- function(fg, ig)
{
	check.svd <- function(m, fg.res)
	{
		mt <- t(as.matrix(m))
		R <- m %*% (mt %*% fg.res$left) - fg.res$left %*% diag(fg.res$values * fg.res$values)
		cat("left-singular vectors:", max(abs(R)), "\n")
		R <- mt %*% (m %*% fg.res$right) - fg.res$right %*% diag(fg.res$values * fg.res$values)
		cat("right-singular vectors:", max(abs(R)), "\n")
		ig.res <- irlba(m, nu=5, nv=5, tol = 1e-12)
		cat("diff on singular values:", sum(abs(fg.res$values) - abs(ig.res$d)), "\n")
	}

	# Get the largest connected component.
	fg.cc <- fg.clusters(fg, mode="weak")
	res <- as.data.frame(table(fg.cc))
	lcc.id <- as.integer(levels(res$fg.cc)[which.max(res$Freq)])
	fg <- fg.fetch.subgraph(fg, which(fg.cc == lcc.id) - 1)
	stopifnot(sum(fg.degree(fg) == 0) == 0)
	ig.matrix <- get.adjacency(ig)[fg.cc == lcc.id, fg.cc == lcc.id]

	# test SVD
	print("test SVD")
	fg.res <- fg.SVD(fg, which="LM", nev=5, ncv=10)
	check.svd(ig.matrix, fg.res)

	# run ASE on the largest connected component

	print("test ASE (adjacency matrix)")
	fg.res <- fg.ASE.igraph(fg, 5, which="A", which.eigen="LM", tol=1.0e-12)
	ig.matrix.tmp <- ig.matrix
	check.svd(ig.matrix, fg.res)

	print("test ASE (AcD)")
	fg.res <- fg.ASE.igraph(fg, 5, which="AcD", which.eigen="LM", c=1/fg.vcount(fg), tol=1.0e-12)
	ig.matrix.tmp <- ig.matrix + diag(fg.degree(fg)) / fg.vcount(fg)
	check.svd(ig.matrix.tmp, fg.res)

	print("test ASE (Laplacian)")
	fg.res <- fg.ASE.igraph(fg, 5, which="L", which.eigen="LM", tol=1.0e-12)
	ig.matrix.tmp <- diag(fg.degree(fg)) - ig.matrix
	check.svd(ig.matrix.tmp, fg.res)

	print("test ASE (normalized Laplacian)")
	fg.res <- fg.ASE.igraph(fg, 5, which="nL", which.eigen="LM", tol=1.0e-12)
	d.matrix <- diag(1 / sqrt(fg.degree(fg)))
	ig.matrix.tmp <- d.matrix %*% (diag(fg.degree(fg)) - ig.matrix) %*% d.matrix
	check.svd(ig.matrix.tmp, fg.res)

	print("test ASE (regularized Laplacian)")
	fg.res <- fg.ASE.igraph(fg, 5, which="nL_tau", which.eigen="LM", tol=1.0e-12, tau=1)
	d.matrix <- diag(1 / sqrt(fg.degree(fg) + 1))
	ig.matrix.tmp <- d.matrix %*% (diag(fg.degree(fg)) - ig.matrix) %*% d.matrix
	check.svd(ig.matrix.tmp, fg.res)
}

test.eigen.undirected <- function(fg, ig)
{
	check.eigen <- function(m, fg.res)
	{
		R <- m %*% fg.res$vectors - fg.res$vectors %*% diag(fg.res$values)
		cat("eigenvectors:", max(abs(R)), "\n")

		multiply <- function(x, extra) {
			(m %*% x)[,1]
		}
		ig.res <- arpack(multiply, sym=TRUE, options=list(n=vcount(ig), nev=5, ncv=10, which="LM"))
		cat("diff on eigen values:", sum(abs(fg.res$values) - abs(ig.res$values)), "\n")
	}
	ig.matrix <- get.adjacency(ig)

	# test eigen
	print("test eigen")
	fg.res <- fg.eigen(fg, which="LM", nev=5, ncv=10)
	check.eigen(ig.matrix, fg.res)

	#test ASE
	print("test ASE (adjacency matrix)")
	fg.res <- fg.ASE.igraph(fg, 5, which="A", which.eigen="LM", tol=1.0e-12)
	check.eigen(ig.matrix, fg.res)

	print("test ASE (AcD)")
	fg.res <- fg.ASE.igraph(fg, 5, which="AcD", which.eigen="LM", c=1/fg.vcount(fg), tol=1.0e-12)
	ig.matrix.tmp <- ig.matrix + diag(fg.degree(fg)) / fg.vcount(fg)
	check.eigen(ig.matrix.tmp, fg.res)

	print("test ASE (Laplacian)")
	fg.res <- fg.ASE.igraph(fg, 5, which="L", which.eigen="LM", tol=1.0e-12)
	ig.matrix.tmp <- diag(fg.degree(fg)) - ig.matrix
	check.eigen(ig.matrix.tmp, fg.res)

	print("test ASE (normalized Laplacian)")
	fg.res <- fg.ASE.igraph(fg, 5, which="nL", which.eigen="LM", tol=1.0e-12)
	d.matrix <- diag(1 / sqrt(fg.degree(fg)))
	ig.matrix.tmp <- d.matrix %*% (diag(fg.degree(fg)) - ig.matrix) %*% d.matrix
	check.eigen(ig.matrix.tmp, fg.res)

	print("test ASE (regularized Laplacian)")
	fg.res <- fg.ASE.igraph(fg, 5, which="nL_tau", which.eigen="LM", tol=1.0e-12, tau=1)
	d.matrix <- diag(1 / sqrt(fg.degree(fg) + 1))
	ig.matrix.tmp <- d.matrix %*% (diag(fg.degree(fg)) - ig.matrix) %*% d.matrix
	check.eigen(ig.matrix.tmp, fg.res)
}

# Kmeans
source("verify.kmeans.R")
test.kmeans(2, 500, 20)
test.with.iris()

# Betweeness
test.betweenness <- function(fg, ig)
{
	fg.btw <- fg.betweenness(fg)
	ig.btw <- betweenness(ig)
	stopifnot(all.equal(ig.btw, fg.btw, tolerance=.001)) # Need tolerance for floats
}

test.weighted <- function(fg, ig)
{
	fg.out <- fg.multiply(fg, rep.int(1, fg.vcount(fg)))
	adj <- get.adjacency(ig, attr='V3', type="both")
	ig.out <- adj %*% rep.int(1, vcount(ig))
	stopifnot(sum(fg.out) == sum(ig.out))
}

# Test on a directed graph.
ig <- read.graph("wiki-Vote1.txt")

fg <- fg.load.graph("wiki-Vote.adj-v4", index="wiki-Vote.index-v4", graph.name="wiki")
test.eigen.directed(fg, ig)

print("run in the standalone mode")
print("load a graph in adjacency list")
fg <- fg.load.graph("wiki-Vote.adj-v4", index="wiki-Vote.index-v4", graph.name="wiki")
test.directed(fg, ig)

cat("\n\n\n")
print("load a graph in edge lists")
fg <- fg.load.graph("wiki-Vote1.txt", graph.name="wiki")
test.directed(fg, ig)
fg.list.graphs()

cat("\n\n\n")
print("load a graph from igraph")
fg <- fg.load.igraph(ig, graph.name="wiki")
test.directed(fg, ig)
fg.list.graphs()

#cat("\n\n\n")
print("run in the SAFS mode")
fg.set.conf("run_test.txt")
fg <- fg.get.graph("wiki-Vote")
test.directed(fg, ig)
fg.list.graphs()

# Now test on an undirected graph
ig <- read.graph("facebook_combined1.txt", directed=FALSE)

fg <- fg.load.graph("facebook.adj-v4", index="facebook.index-v4", graph.name="facebook")
test.eigen.undirected(fg, ig)

print("load a graph in adjacency list")
fg <- fg.load.graph("facebook.adj-v4", index="facebook.index-v4", graph.name="facebook")
test.undirected(fg, ig)

cat("\n\n\n")
print("load a graph in edge lists")
fg <- fg.load.graph("facebook_combined1.txt", directed=FALSE, graph.name="facebook")
test.undirected(fg, ig)
fg.list.graphs()

cat("\n\n\n")
print("load a graph from igraph")
fg <- fg.load.igraph(ig, graph.name="facebook")
test.undirected(fg, ig)
fg.list.graphs()

cat("\n\n\n")
fg <- fg.get.graph("facebook")
test.undirected(fg, ig)
fg.list.graphs()

# Now test on a weighted undirected graph
print("load a weighted graph")
fg <- fg.load.graph("fb-weighted.adj-v4", index="fb-weighted.index-v4", graph.name="fb-weighted")
edges <- read.table("fb-weighted.txt")
ig <- graph.data.frame(edges, directed=FALSE)
test.weighted(fg, ig)
