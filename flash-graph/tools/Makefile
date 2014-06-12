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

OMP_FLAG = -fopenmp
LDFLAGS := -L../../libsafs -lsafs -L.. -lgraph -lrt $(OMP_FLAG) $(LDFLAGS)
CXXFLAGS += -I../../include -I.. -I. $(OMP_FLAG)

all: el2al print_ts_graph graph_merger rmat-gen

el2al: el2al.o ext_mem_vertex_iterator.o ../libgraph.a
	$(CXX) -o el2al el2al.o ext_mem_vertex_iterator.o $(LDFLAGS)

graph_merger: graph_merger.o ext_mem_vertex_iterator.o ../libgraph.a
	$(CXX) -o graph_merger graph_merger.o ext_mem_vertex_iterator.o $(LDFLAGS)

print_ts_graph: print_ts_graph.o ../libgraph.a
	$(CXX) -o print_ts_graph print_ts_graph.o $(LDFLAGS)

rmat-gen: rmat-gen.o
	$(CXX) -o rmat-gen rmat-gen.o

clean:
	rm -f *.d
	rm -f *.o
	rm -f *~
	rm -f el2al
	rm -f graph_merger
	rm -f print_ts_graph
	rm -f rmat-gen

-include $(DEPS) 
