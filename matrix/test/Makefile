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
CXXFLAGS += -I../../include -I../../libcommon -I../../flash-graph -I.. $(OMP_FLAG)
LDFLAGS := $(OMP_FLAG) -L.. -lFMatrix -L../../flash-graph -lgraph -L../../libsafs -lsafs -L../../libcommon -lcommon $(LDFLAGS)
LDFLAGS += -lstxxl -lz -lnuma -laio

all: test-fg_SpMV test-2d_multiply test-mem_matrix test-EM_matrix test-mem_vector	\
	el2al al22d test-sort

test-fg_SpMV: test-fg_SpMV.o ../libFMatrix.a
	$(CXX) -o test-fg_SpMV test-fg_SpMV.o $(LDFLAGS)

test-2d_multiply: test-2d_multiply.o ../libFMatrix.a
	$(CXX) -o test-2d_multiply test-2d_multiply.o $(LDFLAGS)

test-mem_vector: test-mem_vector.o ../libFMatrix.a
	$(CXX) -o test-mem_vector test-mem_vector.o $(LDFLAGS)

test-mem_matrix: test-mem_matrix.o ../libFMatrix.a
	$(CXX) -o test-mem_matrix test-mem_matrix.o $(LDFLAGS)

test-EM_matrix: test-EM_matrix.o ../libFMatrix.a
	$(CXX) -o test-EM_matrix test-EM_matrix.o $(LDFLAGS)

el2al: el2al.o ../libFMatrix.a
	$(CXX) -o el2al el2al.o $(LDFLAGS)

al22d: al22d.o ../libFMatrix.a
	$(CXX) -o al22d al22d.o $(LDFLAGS)

test-sort: test-sort.o ../libFMatrix.a
	$(CXX) -o test-sort test-sort.o $(LDFLAGS)

clean:
	rm -f *.d
	rm -f *.o
	rm -f *~
	rm test-fg_SpMV
	rm test-2d_multiply
	rm test-mem_matrix
	rm test-mem_vector
	rm test-EM_matrix
	rm el2al
	rm al22d
	rm test-sort

-include $(DEPS) 
