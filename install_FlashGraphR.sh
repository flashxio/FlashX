#!/bin/sh

R CMD build flash-graph/Rpkg
if [ -d ~/R ]; then
	path=`find ~/R -name Rcpp.h`
fi
if [ -z "$path" ]; then
	path=`find /usr/local -name Rcpp.h`
fi

matrix_path=`find build -name libmatrix.a`
if [ -n "$matrix_path" ]; then
	echo "find libmatrix.a"
	export FG_EIGEN=1
fi
if [ -n "$path" ]; then
	path=`dirname $path`
	echo $path
	RCPP_INCLUDE=$path FG_DIR=`pwd` FG_LIB=`pwd`/build R CMD INSTALL FlashGraphR_0.2-0.tar.gz
else
	FG_DIR=`pwd` FG_LIB=`pwd`/build R CMD INSTALL FlashGraphR_0.2-0.tar.gz
fi
rm FlashGraphR_0.2-0.tar.gz
