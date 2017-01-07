#!/bin/sh

FlashXDir=$1
FlashGraphRDir=$2

cd $FlashXDir; make clean; cd -
mkdir -p $FlashGraphRDir
rm -R $FlashGraphRDir/*
cp -R $FlashXDir/FlashGraphR/* $FlashGraphRDir
cp -Rf $FlashXDir/flash-graph $FlashGraphRDir/src
mkdir -p $FlashGraphRDir/src/libsafs
mkdir -p $FlashGraphRDir/src/matrix
cp $FlashXDir/libsafs/*.h $FlashGraphRDir/src/libsafs
cp $FlashXDir/matrix/*.h $FlashGraphRDir/src/matrix
cd $FlashGraphRDir/; autoconf; cd -
