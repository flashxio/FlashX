# Copyright 2014 Open Connectome Project (http;//openconnecto.me)
# Written by Da Zheng (zhengda1936@gmail.com)
#
# This file is part of FlashGraph.
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
LDFLAGS := -L.. -lgraph -L../../libsafs -lsafs -L../../libcommon -lcommon -lrt $(OMP_FLAG) -lz $(LDFLAGS)
CXXFLAGS = -I.. -I../../include -I../../libcommon -g -std=c++0x

SOURCE := $(wildcard *.c) $(wildcard *.cpp)
OBJS := $(patsubst %.c,%.o,$(patsubst %.cpp,%.o,$(SOURCE)))
DEPS := $(patsubst %.o,%.d,$(OBJS))

UNITTEST = test-bitmap test-partitioner test-vertex_index

all: $(UNITTEST)

test-bitmap: test-bitmap.o ../libgraph.a
	$(CXX) -o test-bitmap test-bitmap.o $(LDFLAGS)

test-partitioner: test-partitioner.o ../libgraph.a
	$(CXX) -o test-partitioner test-partitioner.o $(LDFLAGS)

test-vertex_index: test-vertex_index.o ../libgraph.a
	$(CXX) -o test-vertex_index test-vertex_index.o $(LDFLAGS)

clean:
	rm -f *.o
	rm -f *.d
	rm -f *~
	rm -f $(UNITTEST)

-include $(DEPS) 
