# Copyright 2013 Da Zheng
#
# This file is part of SA-GraphLib.
#
# SA-GraphLib is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# SA-GraphLib is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with SA-GraphLib.  If not, see <http://www.gnu.org/licenses/>.

include ../Makefile.common

OMP_FLAG = -fopenmp
CXXFLAGS += -I../include $(OMP_FLAG)
LDFLAGS += $(OMP_FLAG)

all: libgraph graph-bfs triangle-counting scan-statistics sssp sstsg unit-test tools graph-stat test cc k-core page-rank diameter-estimator eigensolver

libgraph: $(OBJS)
	rm -f libgraph.a
	ar -cvq libgraph.a $(OBJS)

graph-bfs: libgraph
	$(MAKE) -C graph-bfs

triangle-counting: libgraph
	$(MAKE) -C triangle-counting

scan-statistics: libgraph
	$(MAKE) -C scan-statistics

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

cc: libgraph
	$(MAKE) -C connected-components

k-core: libgraph
	$(MAKE) -C k-core

page-rank: libgraph
	$(MAKE) -C page-rank

diameter-estimator: libgraph
	$(MAKE) -C diameter-estimator

eigensolver: libgraph
	$(MAKE) -C eigensolver

clean:
	rm -f *.d
	rm -f *.o
	rm -f *~
	rm -f libgraph.a
	find -name core -delete
	$(MAKE) clean -C tools
	$(MAKE) clean -C graph-bfs
	$(MAKE) clean -C triangle-counting
	$(MAKE) clean -C scan-statistics
	$(MAKE) clean -C unit-test
	$(MAKE) clean -C sssp
	$(MAKE) clean -C test
	$(MAKE) clean -C sstsg
	$(MAKE) clean -C graph-stat
	$(MAKE) clean -C k-core
	$(MAKE) clean -C weakly-connected-components
	$(MAKE) clean -C page-rank
	$(MAKE) clean -C strongly-connected-components
	$(MAKE) clean -C diameter-estimator
	$(MAKE) clean -C eigensolver

-include $(DEPS) 
