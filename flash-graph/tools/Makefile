# Copyright 2014 Open Connectome Project (http://openconnecto.me)
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

CC = gcc
CXX = g++
OMP_FLAG = -fopenmp
LDFLAGS := -L.. -lgraph -L../../libsafs -lsafs -L../../libcommon -lcommon -lrt -lstxxl $(OMP_FLAG) $(LDFLAGS) -lz
CXXFLAGS += -I../../include -I../../libcommon -I.. -I. $(OMP_FLAG)

all: el2al rmat-gen graph-stat

el2al: el2al.o ../libgraph.a
	$(CXX) -o el2al el2al.o $(LDFLAGS)

print_ts_graph: print_ts_graph.o ../libgraph.a
	$(CXX) -o print_ts_graph print_ts_graph.o $(LDFLAGS)

rmat-gen: rmat-gen.o
	$(CXX) -o rmat-gen rmat-gen.o $(LDFLAGS)

graph-stat: graph-stat.o ../libgraph.a
	$(CXX) -o graph-stat graph-stat.o $(LDFLAGS)

clean:
	rm -f *.d
	rm -f *.o
	rm -f *~
	rm -f el2al
	rm -f print_ts_graph
	rm -f rmat-gen
	rm -f graph-stat

-include $(DEPS) 
