# Copyright 2014 Open Connectome Project (http://openconnecto.me)
# Written by Disa Mhembere (disa@jh.edu)
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

include ../../../Makefile.common

OMP_FLAG = -fopenmp
CXXFLAGS := -g -std=c++0x -I../../../libsafs -I../../../flash-graph -I.. $(OMP_FLAG)
LDFLAGS := $(OMP_FLAG) -L../../eigensolver -leigen -L../../ -lFMatrix -L../../../flash-graph -lgraph -L../../../libsafs -lsafs $(LDFLAGS)
LDFLAGS += -lnuma -lz -laio -lcblas

all: convert_matrix create_matrix


convert_matrix: convert_matrix.o ../../libFMatrix.a
	$(CXX) -o convert_matrix convert_matrix.o $(LDFLAGS)

create_matrix: create_matrix.o ../../libFMatrix.a
	$(CXX) -o create_matrix create_matrix.o $(LDFLAGS)

clean:
	rm -f *.d
	rm -f *.o
	rm -f *~
	rm -f convert_matrix
	rm -f create_matrix

-include $(DEPS) 
