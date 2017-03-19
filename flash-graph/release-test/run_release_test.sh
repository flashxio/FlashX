#!/bin/sh

mkdir -p data
rm -fR data/*

wget http://snap.stanford.edu/data/wiki-Vote.txt.gz
gunzip wiki-Vote.txt.gz
rm -f wiki-Vote.adj*
rm -f wiki-Vote.index*
../utils/el2fg ../conf/run_test.txt wiki-Vote.txt wiki-Vote
sed '/^#/d' wiki-Vote.txt > wiki-Vote1.txt
rm wiki-Vote.txt

wget http://snap.stanford.edu/data/facebook_combined.txt.gz
gunzip facebook_combined.txt.gz
rm -f facebook.adj*
rm -f facebook.index*
../utils/el2fg -u ../conf/run_test.txt facebook_combined.txt facebook
sed '/^#/d' facebook_combined.txt > facebook_combined1.txt

perl add_rand_weight.pl < facebook_combined.txt > fb-weighted.txt
../utils/el2fg -u -t D ../conf/run_test.txt fb-weighted.txt fb-weighted

rm -f facebook_combined.txt
rm wiki-Vote.adj*
rm wiki-Vote.index*
rm facebook.adj*
rm facebook.index*
rm wiki-Vote1.txt
rm facebook_combined1.txt
rm fb-weighted*
