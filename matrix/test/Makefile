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
LDFLAGS += -lstxxl -lz -lnuma -laio -lcblas #-lprofiler

all: test-fg_SpMV test-2d_multiply test-dense_matrix test-mem_vector	\
	el2al al22d al2crs test-sort test-eigen test-dgemm test-mkl_eigen

trilinos: test-anasazi_eigen test-tpetra_multiply

mkl: test-mkl_multiply test-mkl_dgemm

test-fg_SpMV: test-fg_SpMV.o ../libFMatrix.a
	$(CXX) -o test-fg_SpMV test-fg_SpMV.o $(LDFLAGS)

test-2d_multiply: test-2d_multiply.o ../libFMatrix.a
	$(CXX) -o test-2d_multiply test-2d_multiply.o $(LDFLAGS)

test-mem_vector: test-mem_vector.o ../libFMatrix.a
	$(CXX) -o test-mem_vector test-mem_vector.o $(LDFLAGS)

test-dense_matrix: test-dense_matrix.o ../libFMatrix.a
	$(CXX) -o test-dense_matrix test-dense_matrix.o $(LDFLAGS)

el2al: el2al.o ../libFMatrix.a
	$(CXX) -o el2al el2al.o $(LDFLAGS)

al22d: al22d.o fg_utils.o ../libFMatrix.a
	$(CXX) -o al22d al22d.o fg_utils.o $(LDFLAGS)

al2crs: al2crs.o fg_utils.o ../libFMatrix.a
	$(CXX) -o al2crs al2crs.o fg_utils.o $(LDFLAGS)

test-sort: test-sort.o ../libFMatrix.a
	$(CXX) -o test-sort test-sort.o $(LDFLAGS)

LDTRILINOS = -ltpetra -ltpetrakernels -lteuchoscomm -lteuchoskokkoscomm -lteuchosnumerics -lteuchosparameterlist		\
		-lteuchoscore -lanasazi -lkokkoscore -lkokkoscontainers -lteuchoskokkoscompat

test-tpetra_multiply.o: test-tpetra_multiply.cpp
	mpic++ -c test-tpetra_multiply.cpp $(CXXFLAGS)

test-tpetra_multiply: test-tpetra_multiply.o ../libFMatrix.a
	mpic++ -o test-tpetra_multiply test-tpetra_multiply.o $(LDFLAGS) $(LDTRILINOS)

test-mkl_eigen.o: test-mkl_eigen.cpp
	mpic++ -c test-mkl_eigen.cpp $(CXXFLAGS) -I/home/zhengda/intel/mkl/include/

test-mkl_eigen: test-mkl_eigen.o
	mpic++ test-mkl_eigen.o -o test-mkl_eigen $(LDTRILINOS) $(LDFLAGS) -L/home/zhengda/intel/mkl/lib/intel64/ -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lm

test-anasazi_eigen.o: test-anasazi_eigen.cpp
	mpic++ -c test-anasazi_eigen.cpp $(CXXFLAGS)

test-anasazi_eigen: test-anasazi_eigen.o
	mpic++ test-anasazi_eigen.o -o test-anasazi_eigen $(LDTRILINOS) $(LDFLAGS)

test-mkl_multiply.o: test-mkl_multiply.cpp crs_header.h
	$(CXX) -O3 -c test-mkl_multiply.cpp -I/home/zhengda/intel/mkl/include/ -I. -fopenmp

test-mkl_multiply: test-mkl_multiply.o
	$(CXX) -o test-mkl_multiply test-mkl_multiply.o -L/home/zhengda/intel/mkl/lib/intel64/ -lmkl_intel_ilp64 -lmkl_gnu_thread -lmkl_core -lm -fopenmp

test-mkl_dgemm.o: test-dgemm.cpp
	$(CXX) -O3 -c test-dgemm.cpp -o test-mkl_dgemm.o -I/home/zhengda/intel/mkl/include/ -I. -fopenmp -DUSE_MKL

test-mkl_dgemm: test-mkl_dgemm.o
	$(CXX) -o test-mkl_dgemm test-mkl_dgemm.o -L/home/zhengda/intel/mkl/lib/intel64/ -lmkl_intel_ilp64 -lmkl_gnu_thread -lmkl_core -lm -fopenmp -lm

test-dgemm.o: test-dgemm.cpp
	$(CXX) -O3 -c test-dgemm.cpp -I. -fopenmp

test-dgemm: test-dgemm.o
	$(CXX) -o test-dgemm test-dgemm.o -lcblas -lm -fopenmp -lm

test-oblas_dgemm.o: test-dgemm.cpp
	$(CXX) -O3 -c test-dgemm.cpp -o test-oblas_dgemm.o -I. -I/home/zhengda/local/include -fopenmp

test-oblas_dgemm: test-oblas_dgemm.o
	$(CXX) -o test-oblas_dgemm test-oblas_dgemm.o -lopenblas -lm -fopenmp -lm

TRILINOS_LDFLAGS := -lteuchoscomm -lteuchosnumerics -lteuchosparameterlist -lteuchoscore -lanasazi -lanasaziepetra -lepetra
CXXFLAGS += -I../eigensolver

test-eigen: test-eigen.o ../libFMatrix.a ../eigensolver/libeigen.a
	$(CXX) test-eigen.o -o test-eigen -L../eigensolver -leigen $(TRILINOS_LDFLAGS) $(LDFLAGS)

clean:
	rm -f *.d
	rm -f *.o
	rm -f *~
	rm -f test-fg_SpMV
	rm -f test-2d_multiply
	rm -f test-dense_matrix
	rm -f test-mem_vector
	rm -f el2al
	rm -f al22d
	rm -f al2crs
	rm -f test-sort
	rm -f test-tpetra_multiply
	rm -f test-mkl_multiply
	rm -f test-eigen
	rm -f BlockDavidsonTpetra

-include $(DEPS) 
