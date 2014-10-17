require(igraph)
library("FlashGraph")

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

check.vectors <- function(name, v1, v2)
{
	cmp.res <- v1 == v2
	if (sum(cmp.res) != length(cmp.res)) {
		cat(name, " fails the test\n");
	}
}

test.directed <- function(fg, ig)
{
	# test ccoreness
	print("test coreness")
	fg.res <- fg.coreness(fg)
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

	# test degree
	print("test directed degree")
	fg.res <- fg.degree(fg)
	ig.res <- degree(ig)
	check.vectors("degree_test", fg.res, ig.res)

	# test PageRank
	print("test PageRank")
	fg.res <- fg.page.rank(fg)
	ig.res <- page.rank.old(ig, eps=0.01, old=TRUE)
	sum((abs(fg.res - ig.res) / abs(fg.res)) < 0.02)

	# test locality scan
	print("test locality statistics")
	fg.res <- fg.local.scan(fg)
	ig.res <- sapply(graph.neighborhood(ig, 1, mode="all"), ecount)
	check.vectors("local-scan_test", fg.res, ig.res)

	print("test topK locality statistics")
	fg.res <- fg.topK.scan(fg, K=10)
	ig.res <- sort(ig.res, decreasing=TRUE)[1:10]
	check.vectors("topK-scan_test", fg.res$scan, ig.res)
}

test.undirected <- function(fg, ig)
{
	# test triangles
	print("test triangle counting on an undirected graph")
	fg.res <- fg.undirected.triangles(fg)
	ig.res <- adjacent.triangles(ig)
	check.vectors("undirected-triangle_test", fg.res, ig.res)

	# test degree
	print("test undirected degree")
	fg.res <- fg.degree(fg)
	ig.res <- degree(ig)
	check.vectors("degree_test", fg.res, ig.res)
}

fg.set.conf("run_test.txt")
fg <- fg.get.graph("wiki-Vote")
ig <- read.graph("wiki-Vote.txt")
test.directed(fg, ig)

# Now test on an undirected graph
fg <- fg.get.graph("facebook")
ig <- read.graph("facebook_combined.txt", directed=FALSE)
test.undirected(fg, ig)
