# Copyright 2014 Open Connectome Project (http://openconnecto.me)
# Written by Da Zheng (zhengda1936@gmail.com)
#
# This file is part of FlashMatrix.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

include ../../Makefile.common

OMP_FLAG = -fopenmp
CXXFLAGS := -g -std=c++0x -I../../include -I../../libcommon -I../../flash-graph -I.. $(OMP_FLAG)
LDFLAGS := $(OMP_FLAG) -L../ -lFMatrix -L../../flash-graph -lgraph -L../../libsafs -lsafs -L../../libcommon -lcommon $(LDFLAGS)
LDFLAGS += -lstxxl -lnuma -lz -laio

all: test_io_gen test-mem_matrix test-EM_vector test-EM_matrix	\
	test-mem_vector test-data_frame test-mem_vector_vector test-sparse_matrix	\
	test-hilbert test-NUMA_vector test-sorter test-NUMA_dense_matrix

test_io_gen: test_io_gen.o ../libFMatrix.a
	$(CXX) -o test_io_gen test_io_gen.o $(LDFLAGS)

test-mem_matrix: test-mem_matrix.o ../libFMatrix.a
	$(CXX) -o test-mem_matrix test-mem_matrix.o $(LDFLAGS)

test-NUMA_dense_matrix: test-NUMA_dense_matrix.o ../libFMatrix.a
	$(CXX) -o test-NUMA_dense_matrix test-NUMA_dense_matrix.o $(LDFLAGS)

test-mem_vector: test-mem_vector.o ../libFMatrix.a
	$(CXX) -o test-mem_vector test-mem_vector.o $(LDFLAGS)

test-NUMA_vector: test-NUMA_vector.o ../libFMatrix.a
	$(CXX) -o test-NUMA_vector test-NUMA_vector.o $(LDFLAGS)

test-mem_vector_vector: test-mem_vector_vector.o ../libFMatrix.a
	$(CXX) -o test-mem_vector_vector test-mem_vector_vector.o $(LDFLAGS)

test-data_frame: test-data_frame.o ../libFMatrix.a
	$(CXX) -o test-data_frame test-data_frame.o $(LDFLAGS)

test-EM_vector: test-EM_vector.o ../libFMatrix.a
	$(CXX) -o test-EM_vector test-EM_vector.o $(LDFLAGS)

test-EM_matrix: test-EM_matrix.o ../libFMatrix.a
	$(CXX) -o test-EM_matrix test-EM_matrix.o $(LDFLAGS)

test-sparse_matrix: test-sparse_matrix.o ../libFMatrix.a
	$(CXX) -o test-sparse_matrix test-sparse_matrix.o $(LDFLAGS)

test-hilbert: test-hilbert.o ../libFMatrix.a
	$(CXX) -o test-hilbert test-hilbert.o $(LDFLAGS)

test-sorter: test-sorter.o ../libFMatrix.a
	$(CXX) -o test-sorter test-sorter.o $(LDFLAGS)

clean:
	rm -f *.d
	rm -f *.o
	rm -f *~
	rm -f test_io_gen
	rm -f test-mem_matrix
	rm -f test-NUMA_dense_matrix
	rm -f test-mem_vector
	rm -f test-NUMA_vector
	rm -f test-mem_vector_vector
	rm -f test-data_frame
	rm -f test-EM_vector
	rm -f test-EM_matrix
	rm -f test-sparse_matrix
	rm -f test-hilbert
	rm -f test-sorter

-include $(DEPS) 
