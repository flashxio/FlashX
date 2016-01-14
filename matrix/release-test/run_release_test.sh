#!/bin/sh

if [ "$#" -lt 1 ]
then
	echo "test conf_file"
	exit 1
fi

conf_file=$1

echo "unit test on data frame"
../unit-test/test-data_frame $conf_file
echo "unit test on external-memory matrix"
../unit-test/test-EM_matrix $conf_file
echo "unit test on external-memory vector"
../unit-test/test-EM_vector $conf_file
echo "unit test on external-memory vector vector"
../unit-test/test-EM_vector_vector $conf_file
echo "unit test on local matrix store"
../unit-test/test-local_matrix_store
echo "unit test on in-memory matrix store"
../unit-test/test-mem_matrix_store
echo "unit test on in-memory special matrix store"
../unit-test/test-special_matrix_store
echo "unit test on in-memory vector"
../unit-test/test-mem_vector
echo "unit test on in-memory vector of vectors"
../unit-test/test-mem_vector_vector
echo "unit test on sparse matrix"
../unit-test/test-sparse_matrix
echo "unit test on NUMA vector"
../unit-test/test-NUMA_vector $conf_file
echo "unit test on sorter"
../unit-test/test-sorter
echo "unit test on NUMA dense matrix"
../unit-test/test-NUMA_dense_matrix $conf_file
echo "unit test on dense matrix"
../unit-test/test-dense_matrix $conf_file

rm -f wiki-Vote.txt.gz
wget http://snap.stanford.edu/data/wiki-Vote.txt.gz
gunzip wiki-Vote.txt.gz
sed '/^#/d' wiki-Vote.txt > wiki-Vote1.txt
rm wiki-Vote.txt

R --no-save < verify.R
#R -d valgrind --vanilla < verify.R
#R -d gdb --vanilla

rm wiki-Vote1.txt
