require(igraph)
require(argparse)

parser <- ArgumentParser(description="Test igraph wiki result vs FlashGraph")
parser$add_argument("vfn", help="The FG vector written to file using FG_vector<T>::ptr->to_file()")
parser$add_argument("-w", "--wikifn", default="/mnt/nfs2/graphs/wiki-Vote.txt", help="The file location of the wiki-Vote edgelist file on disk")

result <- parser$parse_args()

g <- read.graph(result$wikifn)
ig.coreness <- graph.coreness(g, mode="all")
suppressWarnings(fg.coreness <- as.numeric(read.delim(result$vfn, header = FALSE, sep = " ")[1,])[1])
fg.coreness <- fg.coreness[1:length(fg.coreness)-1]

if (all(ig.coreness == fg.coreness)) 
	cat("\nSUCCESS!! igraph and flashgraph coreness values match!\n") else
	cat("\nMISMATCH! igraph and flashgraph values DO NOT match!\n")

