num.clusters <- 1000
num.eigens <- 40
library(FlashR)
fg.set.conf("flash-graph/conf/run_ssds.txt")
vecs <- fm.read.obj(file="/mnt/store_raid/zhengda/friendster-lcc.acd.LM40.1e-6.evecs")
col3 <- fm.get.cols(vecs, 1:3)
r.col3 <- fm.conv.FM2R(col3)
vals <- fm.conv.FM2R(fm.read.obj(file="/mnt/store_raid/zhengda/friendster-lcc.acd.LM40.1e-6.evals"))[1:3]
ase <- r.col3 %*% diag(abs(as.vector(vals)))
sXhat <- ase / sqrt(rowSums(ase^2))
kmeans.res <- fg.kmeans(sXhat, num.clusters, max.iters=100, init="random")
fg <- fg.load.graph("/mnt/store_raid/zhengda/friendster-lcc.adj-v4",
"/mnt/store_raid/zhengda/friendster-lcc.index-v4")

large.ccs <- c()
clusters <- kmeans.res$clusters
for (cluster.id in 1:num.clusters) {
	clusterV <- which(kmeans.res$cluster==cluster.id)
	if (length(clusterV) == 0)
		next
	sub.fg <- fg.fetch.subgraph(fg, vertices=clusterV - 1, compress=FALSE)
	cat("subg #V:", fg.vcount(sub.fg), ", #E:", fg.ecount(sub.fg), "\n")
	sub.cc <- fg.clusters(sub.fg, mode="weak")
	clusters[clusterV] <- sub.cc[clusterV]

	counts <- as.data.frame(table(sub.cc[sub.cc >= 0]))
	lcc.id <- as.integer(levels(counts$Var1)[which.max(counts$Freq)])
	if (max(counts$Freq) > 100)
		large.ccs <- c(large.ccs, lcc.id)
	cat("subg #components:", length(counts$Var1), ", lcc size:",
		sum(sub.cc == lcc.id), ", #isolates:", sum(sub.cc == -1), "\n")
	# lcc.fg <- fg.fetch.subgraph(sub.fg, vertices=)
	cat("There are ", length(unique(clusters)), " clusters\n")
}

# This can be used to collect the cluster membership from the vertices
# in the components.
tmp.clusters <- rep.int(-1, fg.vcount(fg))
max.cid <- 0
for (large.cc in large.ccs) {
	clusterV <- which(clusters==large.cc)
	clusters[clusterV] <- -1
	sub.fg <- fg.fetch.subgraph(fg, vertices=clusterV - 1)
	# The cluster IDs are only valid in the subgraph.
	sub.clusters <- spec.cluster(sub.fg)
	# Now the cluster IDs are globally unique.
	tmp.clusters[clusterV] <- sub.clusters + max.cid
	max.cid <- max.cid + fg.vcount(sub.fg)
}
clusters[clusters != -1] <- clusters[clusters != -1] + max.cid
clusters[tmp.clusters != -1] <- tmp.clusters[tmp.clusters != -1]

spec.cluster <- function(fg)
{
	eigen <- fg.ASE.igraph(fg, num.eigen=num.eigens, which="AcD", which.eigen="LM")
	ase <- eigen$vectors %*% diag(abs(as.vector(eigen$values)))
	sXhat <- ase / sqrt(rowSums(ase^2))
	kmeans.res <- fg.kmeans(sXhat, num.clusters, max.iters=100, init="random")

	clusters <- kmeans.res$clusters
	for (cluster.id in 1:num.clusters) {
		clusterV <- which(kmeans.res$cluster==cluster.id)
		if (length(clusterV) == 0)
			next
		sub.fg <- fg.fetch.subgraph(fg, vertices=clusterV - 1, compress=FALSE)
		cat("subg #V:", fg.vcount(sub.fg), ", #E:", fg.ecount(sub.fg), "\n")
		sub.cc <- fg.clusters(sub.fg, mode="weak")
		clusters[clusterV] <- sub.cc[clusterV]
	}
	clusters
}
