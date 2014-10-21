library("FlashGraph")
library("igraph")

graph.name <- "wiki-Vote"
cc.type <- "weak"

fg <- fg.get.graph(graph.name)
cc <- fg.clusters(fg, mode=cc.type)
cc.counts <- as.data.frame(table(cc))
max.cc <- max(cc.counts$Freq)
cat("There are", nrow(cc.counts), "components\n")
cat("The largest component has", max.cc, "vertices\n")

fetch.cc <- function(cid)
{
	# Vertex ID starts with 0.
	vids <- which(cc == cid) - 1
	fg.fetch.subgraph(fg, vids)
}

# fetch each component from the graph but the largest one.
stat.cc <- function(x)
{
	cc.id <- x[1]
	cc.size <- x[2]
	if (cc.id != -1 && cc.size != max.cc) {
		subg <- fetch.cc(cc.id)
		print(sprintf("component %d has %d vertices and %d edges\n",
					  as.integer(cc.id), vcount(subg), ecount(subg)))
	}
}

apply(cc.counts, 1, stat.cc)
