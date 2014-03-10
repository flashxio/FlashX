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

include ../../Makefile.common

OMP_FLAG = -fopenmp
LDFLAGS := -L../../libsafs -lsafs -L.. -lgraph -lrt $(OMP_FLAG) $(LDFLAGS)
CXXFLAGS += -I../../include -I.. -I. $(OMP_FLAG)

all: el2al print_ts_graph graph_merger

el2al: el2al.o graph.o ../libgraph.a
	$(CXX) -o el2al el2al.o graph.o $(LDFLAGS)

graph_merger: graph_merger.o graph.o ../libgraph.a
	$(CXX) -o graph_merger graph_merger.o graph.o $(LDFLAGS)

print_ts_graph: print_ts_graph.o ../libgraph.a
	$(CXX) -o print_ts_graph print_ts_graph.o $(LDFLAGS)

clean:
	rm -f *.d
	rm -f *.o
	rm -f *~
	rm -f el2al
	rm -f graph_merger
	rm -f print_ts_graph

-include $(DEPS) 
