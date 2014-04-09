# Copyright 2013 Disa Mhembere
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
LDFLAGS := -L.. -lgraph -L../../libsafs -lsafs -lrt $(OMP_FLAG) $(LDFLAGS)
CXXFLAGS += -I../../include -I.. -I. $(OMP_FLAG)

all: page-rank page-rank2

page-rank: page_rank.o ../libgraph.a
	$(CXX) -o page-rank	page_rank.o $(LDFLAGS)

page-rank2: page_rank2.o ../libgraph.a
	$(CXX) -o page-rank2 page_rank2.o $(LDFLAGS)

clean:
	rm -f *.d
	rm -f *.o
	rm -f *~
	rm -f page-rank

-include $(DEPS) 
