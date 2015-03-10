#!/bin/sh

../unit-test/test-data_frame
../unit-test/test-EM_matrix
../unit-test/test-EM_vector
../unit-test/test-mem_matrix
../unit-test/test-mem_vector
../unit-test/test-mem_vector_vector

wget http://snap.stanford.edu/data/wiki-Vote.txt.gz
gunzip wiki-Vote.txt.gz
sed '/^#/d' wiki-Vote.txt > wiki-Vote1.txt
rm wiki-Vote.txt

R --no-save < verify.R

rm wiki-Vote1.txt
