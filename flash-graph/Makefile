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
CXXFLAGS += -I../include -I../libcommon $(OMP_FLAG)
LDFLAGS += $(OMP_FLAG)

all: libgraph unit-test tools test matrix libgraph-algs test-algs test-matrix

libgraph-algs: libgraph
	$(MAKE) -C libgraph-algs

libgraph: $(OBJS)
	rm -f libgraph.a
	ar -cvq libgraph.a $(OBJS)

tools: libgraph
	$(MAKE) -C tools

unit-test: libgraph
	$(MAKE) -C unit-test

test: libgraph
	$(MAKE) -C test

matrix: libgraph
	$(MAKE) -C matrix

test-algs: libgraph libgraph-algs matrix
	$(MAKE) -C test-algs

LOL: libgraph matrix
	$(MAKE) -C LOL

test-matrix: libgraph matrix
	$(MAKE) -C test-matrix

clean:
	rm -f *.d
	rm -f *.o
	rm -f *~
	rm -f libgraph.a
	find -name core -delete
	$(MAKE) clean -C tools
	$(MAKE) clean -C unit-test
	$(MAKE) clean -C test
	$(MAKE) clean -C matrix
	$(MAKE) clean -C libgraph-algs
	$(MAKE) clean -C test-algs
	$(MAKE) clean -C LOL
	$(MAKE) clean -C test-matrix

-include $(DEPS) 
