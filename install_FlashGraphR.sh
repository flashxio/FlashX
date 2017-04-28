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

R CMD build FlashGraphR
version=`R --version | grep "R version" | awk '{print $3}' | awk -F. '{print $1 "." $2}'`
if [ -d ~/R ]; then
	path=`find ~/R -name Rcpp.h | grep ${version}`
fi
if [ -z "$path" ]; then
	path=`find /usr/local -name Rcpp.h`
fi

echo "search for FlashR"
libpaths=`R -e 'cat(.libPaths(), "\n")' --no-save --slave`
FlashR_dir=""
for x in $libpaths
do
	if [ -f "$x/FlashR/libs/FlashR.so" ]; then
		FlashR_dir="$x/FlashR/libs/"
	fi
done

if [ "$FlashR_dir" = "" ]; then
	echo "cannot find FlashR"
	exit 1
fi

export FlashR_dir
fg_lib=`pwd`/build
if [ -n "$path" ]; then
	path=`dirname $path`
	echo $path
	RCPP_INCLUDE=$path FG_DIR=`pwd` FG_LIB=$fg_lib R CMD INSTALL --no-configure FlashGraphR_0.1-0.tar.gz
else
	FG_DIR=`pwd` FG_LIB=$fg_lib R CMD INSTALL --no-configure FlashGraphR_0.1-0.tar.gz
fi
rm FlashGraphR_0.1-0.tar.gz
