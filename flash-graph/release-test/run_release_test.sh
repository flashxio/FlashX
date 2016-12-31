#!/bin/sh

mkdir -p data
rm -fR data/*

wget http://snap.stanford.edu/data/wiki-Vote.txt.gz
gunzip wiki-Vote.txt.gz
rm wiki-Vote.adj*
rm wiki-Vote.index*
../../matrix/test/el2al ../conf/run_test.txt wiki-Vote.txt wiki-Vote
sed '/^#/d' wiki-Vote.txt > wiki-Vote1.txt
rm wiki-Vote.txt

wget http://snap.stanford.edu/data/facebook_combined.txt.gz
gunzip facebook_combined.txt.gz
rm facebook.adj*
rm facebook.index*
../../matrix/test/el2al -u ../conf/run_test.txt facebook_combined.txt facebook
sed '/^#/d' facebook_combined.txt > facebook_combined1.txt

#perl add_rand_weight.pl < facebook_combined.txt > fb-weighted.txt
#../../matrix/test/el2al -u -t I ../conf/run_test.txt fb-weighted.txt fb-weighted

rm facebook_combined.txt

OMP_NUM_THREADS=1 R --no-save < verify.vs.igraph.R

rm wiki-Vote.adj*
rm wiki-Vote.index*
rm facebook.adj*
rm facebook.index*
rm wiki-Vote1.txt
rm facebook_combined1.txt
rm fb-weighted*

R --no-save < verify.large.R
