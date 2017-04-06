library(testthat)
library(FlashR)

ret <- test_file("Rpkg/inst/tests/test_dense.R")
msgs <- capture.output(ret)
tests <- length(msgs) / 3
fails <- sum(as.integer(sapply(strsplit(msgs[(tests+2):(2*tests)], " +"),
					  function(x) x[length(x)-1])))
q(status=fails)
