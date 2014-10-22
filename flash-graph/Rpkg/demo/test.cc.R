library("FlashGraph")
library("igraph")

test.fg.cc <- function(fg, cc.type)
{
	cc <- fg.clusters(fg, mode=cc.type)
	cc.counts <- as.data.frame(table(cc))
	max.cc <- max(cc.counts$Freq)
	cat("There are", nrow(cc.counts), cc.type, "connected components\n")
	cat("The largest component has", max.cc, "vertices\n")

	fg.fetch.cc <- function(cid)
	{
		# Vertex ID starts with 0.
		vids <- which(cc == cid) - 1
		fg.fetch.subgraph(fg, vids)
	}

	# fetch each component from the graph but the largest one.
	fg.stat.cc <- function(x)
	{
		cc.id <- x[1]
		cc.size <- as.integer(x[2])
		if (cc.id != -1 && cc.size > 1) {
			subg <- fg.fetch.cc(cc.id)
			ecount(subg)
		}
	}
	apply(cc.counts, 1, fg.stat.cc)
}

test.ig.cc <- function(ig, cc.type)
{
	ig.cc <- clusters(ig, mode=cc.type)
#	isolates <- which(ig.cc$csize==1)
#	ig.cc$membership[ig.cc$membership %in% isolates]=-1
	ig.cc.counts <- as.data.frame(table(ig.cc$membership))
	ig.max.cc <- max(ig.cc.counts$Freq)
	cat("There are", nrow(ig.cc.counts), cc.type, "connected components\n")
	cat("The largest component has", ig.max.cc, "vertices\n")

	ig.fetch.cc <- function(cid)
	{
		# Vertex ID starts with 0.
		vids <- which(ig.cc$membership == cid)
		induced.subgraph(ig, vids)
	}

	# fetch each component from the graph but the largest one.
	ig.stat.cc <- function(x)
	{
		cc.id <- x[1]
		cc.size <- as.integer(x[2])
		if (cc.id != -1 && cc.size > 1) {
			subg <- ig.fetch.cc(cc.id)
			ecount(subg)
		}
	}
	apply(ig.cc.counts, 1, ig.stat.cc)
}

cc.type <- "weak"
graph.name <- "wiki-Vote"
fg <- fg.get.graph(graph.name)
ig <- read.graph("wiki-Vote.txt")
fg.res <- test.fg.cc(fg, cc.type)
ig.res <- test.ig.cc(ig, cc.type)

unlist(ig.res) == unlist(fg.res)
