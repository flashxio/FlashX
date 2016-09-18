#!/bin/sh

FlashXDir=$1
FlashRDir=$2

cd $FlashXDir; make clean; cd -
cp -R $FlashXDir/Rpkg $FlashRDir
cp -R $FlashXDir/flash-graph $FlashRDir/src
cp -R $FlashXDir/libsafs $FlashRDir/src
cp -R $FlashXDir/matrix $FlashRDir/src
cd $FlashRDir/; aclocal; autoconf; cd -
