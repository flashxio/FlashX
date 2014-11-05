require(igraph)
library("FlashGraphR")
library("irlba")

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
	stopifnot(sum(cmp.res) == length(cmp.res))
}

test.directed <- function(fg, ig)
{
	# test degree
	# this can be used to test the correctness of the generated graph.
	print("test directed degree")
	time1 <- system.time(fg.res <- fg.degree(fg))
	cat("FG:", time1, "\n")
	time2 <- system.time(ig.res <- degree(ig))
	cat("IG:", time2, "\n")
	check.vectors("degree_test", fg.res, ig.res)

	# test ccoreness
	print("test coreness")
	time1 <- system.time(fg.res <- fg.coreness(fg))
	cat("FG:", time1, "\n")
	time2 <- system.time(ig.res <- graph.coreness(ig, mode="all"))
	cat("IG:", time2, "\n")
	check.vectors("coreness_test", fg.res, ig.res)

	# test WCC
	print("test WCC")
	time1 <- system.time(fg.res <- fg.clusters(fg, mode="weak"))
	cat("FG:", time1, "\n")
	time2 <- system.time(ig.res <- clusters(ig, mode="weak")$membership)
	cat("IG:", time2, "\n")
	verify.cc(fg.res, ig.res)

	# test SCC
	print("test SCC")
	time1 <- system.time(fg.res <- fg.clusters(fg, mode="strong"))
	cat("FG:", time1, "\n")
	time2 <- system.time(ig.res <- clusters(ig, mode="strong")$membership)
	cat("IG:", time2, "\n")
	verify.cc(fg.res, ig.res)

	# test PageRank
	print("test PageRank")
	time1 <- system.time(fg.res <- fg.page.rank(fg))
	cat("FG:", time1, "\n")
	time2 <- system.time(ig.res <- page.rank.old(ig, eps=0.01, old=TRUE))
	cat("IG:", time2, "\n")
	num <- sum((abs(fg.res - ig.res) / abs(fg.res)) < 0.02)
	cat("# vertices whose PR diff <= 2% is", num, ", # vertices:", vcount(ig))

	# test locality scan
	print("test locality statistics")
	time1 <- system.time(fg.res <- fg.local.scan(fg))
	cat("FG:", time1, "\n")
	time2 <- system.time(ig.res <- sapply(graph.neighborhood(ig, 1, mode="all"), ecount))
	cat("IG:", time2, "\n")
	check.vectors("local-scan_test", fg.res, ig.res)

	# test topK scan
	print("test topK locality statistics")
	time1 <- system.time(fg.res <- fg.topK.scan(fg, K=10))
	cat("FG:", time1, "\n")
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

	# test SVD
	print("test SVD")
	fg.res <- fg.SVD(fg, which="LM", nev=5, ncv=10)
	ig.res <- irlba(ig.matrix, nu=5, nv=0)
	cat("diff on singular values:", sum(abs(fg.res$values) - abs(ig.res$d)), "\n")
	cat("diff on left-singular vectors:", sum(abs(fg.res$vectors) - abs(ig.res$u)), "\n")
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
	time1 <- system.time(fg.res <- fg.undirected.triangles(fg))
	cat("FG:", time1, "\n")
	time2 <- system.time(ig.res <- adjacent.triangles(ig))
	cat("IG:", time2, "\n")
	check.vectors("undirected-triangle_test", fg.res, ig.res)

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

	# test eigen
	print("test eigen")
	fg.res <- fg.eigen(fg, which="LM", nev=5, ncv=10)
	multiply <- function(x, extra) {
		(ig.matrix %*% x)[,1]
	}
	ig.res <- arpack(multiply, sym=TRUE, options=list(n=vcount(ig), nev=5, ncv=10, which="LM"))
	cat("diff on eigen values:", sum(abs(fg.res$values) - abs(ig.res$values)), "\n")
	cat("diff on eigen vectors:", sum(abs(fg.res$vectors) - abs(ig.res$vectors)), "\n")
}

# Test on a directed graph.
ig <- read.graph("wiki-Vote1.txt")

print("run in the standalone mode")
print("load a graph in adjacency list")
fg <- fg.load.graph("wiki-Vote.adj-v4", index="wiki-Vote.index-v4")
test.directed(fg, ig)

cat("\n\n\n")
print("load a graph in edge lists")
fg <- fg.load.graph("wiki-Vote1.txt")
test.directed(fg, ig)
fg.list.graphs()

cat("\n\n\n")
print("load a graph from igraph")
fg <- fg.load.graph(ig)
test.directed(fg, ig)
fg.list.graphs()

cat("\n\n\n")
print("run in the SAFS mode")
fg.set.conf("run_test.txt")
fg <- fg.get.graph("wiki-Vote")
test.directed(fg, ig)
fg.list.graphs()

# Now test on an undirected graph
ig <- read.graph("facebook_combined1.txt", directed=FALSE)

print("load a graph in adjacency list")
fg <- fg.load.graph("facebook.adj-v4", index="facebook.index-v4")
test.undirected(fg, ig)

cat("\n\n\n")
print("load a graph in edge lists")
fg <- fg.load.graph("facebook_combined1.txt", directed=FALSE)
test.undirected(fg, ig)
fg.list.graphs()

cat("\n\n\n")
print("load a graph from igraph")
fg <- fg.load.graph(ig)
test.undirected(fg, ig)
fg.list.graphs()

cat("\n\n\n")
fg <- fg.get.graph("facebook")
test.undirected(fg, ig)
fg.list.graphs()
