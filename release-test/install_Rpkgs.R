r <- getOption("repos")
r["CRAN"] <- "http://cran.cnr.berkeley.edu/"
options(repos=r)
install.packages("testthat")
install.packages("Rcpp")
install.packages("RSpectra")
