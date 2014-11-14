library(FlashGraphR)

fg.set.conf("run_test_large.txt")

stopifnot(fg.exist.graph("twitter"))
fg <- fg.get.graph("twitter")

count.cc <- function(cc)
{
	sum(cc != -1)
}

print("test WCC")
system.time(fg.res <- fg.clusters(fg, mode="weak"))
fg.res.counts <- as.data.frame(table(fg.res))
cat("largest WCC has", max(fg.res.counts$Freq), "vertices, there are",
	count.cc(fg.res.counts$fg.res), "components\n")
stopifnot(max(fg.res.counts$Freq) == sum(fg.degree(fg) > 0))
stopifnot(count.cc(fg.res.counts$fg.res) == 1)

print("test SCC")
system.time(fg.res <- fg.clusters(fg, mode="strong"))
fg.res.counts <- as.data.frame(table(fg.res))
cat("largest SCC has", max(fg.res.counts$Freq), "vertices, there are",
	count.cc(fg.res.counts$fg.res), "components\n")
stopifnot(max(fg.res.counts$Freq) == 33479734)
stopifnot(count.cc(fg.res.counts$fg.res) == 8044728)

print("test topK scan")
system.time(fg.res <- fg.topK.scan(fg, K=1))
cat("scan statistics:", fg.res$scan, "\n")
stopifnot(fg.res$vid == 11915432)
stopifnot(fg.res$scan == 114616908)
scan.stat <- fg.res$scan

print("A * x");
x <- rep(1, times=fg.vcount(fg))
system.time(fg.res <- fg.multiply(fg, x))
cat("A * x =", sum(fg.res), "\n")
stopifnot(sum(fg.res) == fg.ecount(fg))

print("t(A) * x");
system.time(fg.res <- fg.multiply(fg, x, TRUE))
cat("t(A) * x =", sum(fg.res), "\n")
stopifnot(sum(fg.res) == fg.ecount(fg))

stopifnot(fg.exist.graph("friendster"))
fg <- fg.get.graph("friendster")

print("test triangle counting on an undirected graph")
system.time(fg.res <- fg.undirected.triangles(fg))
stopifnot(sum(as.numeric(fg.res))/3 == 4173724142)

print("A * x");
x <- rep(1, times=fg.vcount(fg))
system.time(fg.res <- fg.multiply(fg, x))
cat("A * x =", sum(fg.res), "\n")
stopifnot(sum(fg.res) == fg.ecount(fg) * 2)

print("t(A) * x");
system.time(fg.res <- fg.multiply(fg, x, TRUE))
cat("t(A) * x =", sum(fg.res), "\n")
stopifnot(sum(fg.res) == fg.ecount(fg) * 2)
