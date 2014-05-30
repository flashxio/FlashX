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

include ../Makefile.common

OMP_FLAG = -fopenmp
CXXFLAGS += -I../include $(OMP_FLAG)
LDFLAGS += $(OMP_FLAG)

all: libgraph graph-bfs sssp sstsg unit-test tools graph-stat test k-core page-rank diameter-estimator eigensolver matrix libgraph-algs test-algs

libgraph-algs: libgraph
	$(MAKE) -C libgraph-algs

libgraph: $(OBJS)
	rm -f libgraph.a
	ar -cvq libgraph.a $(OBJS)

graph-bfs: libgraph
	$(MAKE) -C graph-bfs

sssp: libgraph
	$(MAKE) -C sssp

sstsg: libgraph
	$(MAKE) -C sstsg

tools: libgraph
	$(MAKE) -C tools

unit-test: libgraph
	$(MAKE) -C unit-test

graph-stat: libgraph
	$(MAKE) -C graph-stat

test: libgraph
	$(MAKE) -C test

k-core: libgraph
	$(MAKE) -C k-core

page-rank: libgraph
	$(MAKE) -C page-rank

diameter-estimator: libgraph
	$(MAKE) -C diameter-estimator

eigensolver: libgraph
	$(MAKE) -C eigensolver

matrix: libgraph
	$(MAKE) -C matrix

test-algs: libgraph
	$(MAKE) -C test-algs

clean:
	rm -f *.d
	rm -f *.o
	rm -f *~
	rm -f libgraph.a
	find -name core -delete
	$(MAKE) clean -C tools
	$(MAKE) clean -C graph-bfs
	$(MAKE) clean -C unit-test
	$(MAKE) clean -C sssp
	$(MAKE) clean -C test
	$(MAKE) clean -C sstsg
	$(MAKE) clean -C graph-stat
	$(MAKE) clean -C k-core
	$(MAKE) clean -C connected-components
	$(MAKE) clean -C page-rank
	$(MAKE) clean -C diameter-estimator
	$(MAKE) clean -C eigensolver
	$(MAKE) clean -C matrix
	$(MAKE) clean -C libgraph-algs
	$(MAKE) clean -C test-algs

-include $(DEPS) 
