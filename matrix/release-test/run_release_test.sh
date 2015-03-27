#!/bin/sh

echo "unit test on data frame"
../unit-test/test-data_frame
echo "unit test on external-memory matrix"
../unit-test/test-EM_matrix
echo "unit test on external-memory vector"
../unit-test/test-EM_vector
echo "unit test on in-memory matrix"
../unit-test/test-mem_matrix
echo "unit test on in-memory vector"
../unit-test/test-mem_vector
echo "unit test on in-memory vector of vectors"
../unit-test/test-mem_vector_vector
echo "unit test on sparse matrix"
../unit-test/test-sparse_matrix
echo "unit test on NUMA vector"
../unit-test/test-NUMA_vector
echo "unit test on sorter"
../unit-test/test-sorter

rm -f wiki-Vote.txt.gz
wget http://snap.stanford.edu/data/wiki-Vote.txt.gz
gunzip wiki-Vote.txt.gz
sed '/^#/d' wiki-Vote.txt > wiki-Vote1.txt
rm wiki-Vote.txt

R --no-save < verify.R

rm wiki-Vote1.txt
