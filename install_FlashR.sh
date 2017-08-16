#!/bin/sh

if [ `numactl -H | grep "nodes" | awk '{print $2}'` -gt 1 ]; then
	export ENABLE_NUMA=1
fi

mkdir -p build
cd build
cmake ..;
if [ $? != 0 ]; then
	exit
fi
make
cd ..

R CMD build FlashR
version=`R --version | grep "R version" | awk '{print $3}' | awk -F. '{print $1 "." $2}'`
if [ -d ~/R ]; then
	path=`find ~/R -name Rcpp.h | grep ${version}`
fi
if [ -z "$path" ]; then
	path=`find /usr/local -name Rcpp.h`
fi

fg_lib=`pwd`/build
eigen_path=`find $fg_lib/matrix -name libeigen.a`
if [ -n "$eigen_path" ]; then
	echo "find libeigen.a"
	export ENABLE_TRILINOS=1
fi
if [ -n "$path" ]; then
	path=`dirname $path`
	echo $path
	RCPP_INCLUDE=$path FG_DIR=`pwd` FG_LIB=$fg_lib R CMD INSTALL --no-configure FlashR_0.1-0.tar.gz
else
	FG_DIR=`pwd` FG_LIB=$fg_lib R CMD INSTALL --no-configure FlashR_0.1-0.tar.gz
fi
ret=$?
rm FlashR_0.1-0.tar.gz
exit $ret
