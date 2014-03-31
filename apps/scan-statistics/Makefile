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

all: topK-scan local-scan overlap

topK-scan: scan_graph.o topK_scan_graph.o ../libgraph.a
	$(CXX) -o topK-scan scan_graph.o topK_scan_graph.o $(LDFLAGS)

local-scan: scan_graph.o local_scan_graph.o ../libgraph.a
	$(CXX) -o local-scan scan_graph.o local_scan_graph.o $(LDFLAGS)

overlap: overlap.o ../libgraph.a
	$(CXX) -o overlap overlap.o $(LDFLAGS)

clean:
	rm -f *.d
	rm -f *.o
	rm -f *~
	rm -f topK-scan
	rm -f local-scan
	rm -f overlap

-include $(DEPS) 
