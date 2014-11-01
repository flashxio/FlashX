#!/bin/sh

R CMD build flash-graph/Rpkg
if [ -d ~/R ]; then
	path=`find ~/R -name Rcpp.h`
fi
if [ -n "$path" ]; then
	path=`dirname $path`
	echo $path
	RCPP_INCLUDE=$path FG_DIR=`pwd` FG_LIB=`pwd`/build R CMD INSTALL FlashGraph_0.1-0.tar.gz
else
	FG_DIR=`pwd` FG_LIB=`pwd`/build R CMD INSTALL FlashGraph_0.1-0.tar.gz
fi
rm FlashGraph_0.1-0.tar.gz
