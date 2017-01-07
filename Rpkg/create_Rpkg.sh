#!/bin/sh

FlashXDir=$1
FlashRDir=$2

cd $FlashXDir; make clean; cd -
mkdir -p $FlashRDir
rm -R $FlashRDir/*
cp -R $FlashXDir/Rpkg/* $FlashRDir
cp -Rf $FlashXDir/libsafs $FlashRDir/src
cp -Rf $FlashXDir/matrix $FlashRDir/src
cd $FlashRDir/; aclocal; autoconf; cd -
