#!/bin/sh

mkdir -p data

wget http://snap.stanford.edu/data/wiki-Vote.txt.gz
gunzip wiki-Vote.txt.gz
rm wiki-Vote.adj*
rm wiki-Vote.index*
../tools/el2al -v -w wiki-Vote.adj wiki-Vote.index wiki-Vote.txt
../../utils/SAFS-util run_test.txt load wiki-Vote-v4 wiki-Vote.adj-v4
../../utils/SAFS-util run_test.txt load wiki-Vote-index-v4 wiki-Vote.index-v4
sed '/^#/d' wiki-Vote.txt > wiki-Vote1.txt
rm wiki-Vote.txt

wget http://snap.stanford.edu/data/facebook_combined.txt.gz
gunzip facebook_combined.txt.gz
rm facebook.adj*
rm facebook.index*
../tools/el2al -v -u -w facebook.adj facebook.index facebook_combined.txt
../../utils/SAFS-util run_test.txt load facebook-v4 facebook.adj-v4
../../utils/SAFS-util run_test.txt load facebook-index-v4 facebook.index-v4
sed '/^#/d' facebook_combined.txt > facebook_combined1.txt
rm facebook_combined.txt

OMP_NUM_THREADS=1 R --no-save < verify.vs.igraph.R

rm wiki-Vote.adj*
rm wiki-Vote.index*
rm facebook.adj*
rm facebook.index*
rm wiki-Vote1.txt
rm facebook_combined1.txt

#R --no-save < verify.large.R
