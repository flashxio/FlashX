# Copyright 2016 Open Connectome Project (http://openconnecto.me)
# Written by Da Zheng (zhengda1936@gmail.com)
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

include ../../Makefile.common

CXXFLAGS += -I../../libsafs -I.. -I../../matrix
LDFLAGS := -L../ -lgraph -L../../matrix -lFMatrix -L../../libsafs -lsafs $(LDFLAGS)
LDFLAGS += -lz -lcblas #-lprofiler

all: el2fg fg2fm fg2crs fg_lcc csr2fg sbm

el2fg: el2fg.o ../libgraph.a
	$(CXX) -o el2fg el2fg.o $(LDFLAGS)

csr2fg: csr2fg.o ../libgraph.a
	$(CXX) -o csr2fg csr2fg.o $(LDFLAGS)

fg2fm: fg2fm.o ../libgraph.a
	$(CXX) -o fg2fm fg2fm.o $(LDFLAGS)

fg2crs: fg2crs.o ../libgraph.a
	$(CXX) -o fg2crs fg2crs.o $(LDFLAGS)

fg_lcc: fg_lcc.o ../libgraph.a
	$(CXX) -o fg_lcc fg_lcc.o -L../../flash-graph/libgraph-algs -lgraph-algs $(LDFLAGS)

sbm: sbm.o ../libgraph.a
	$(CXX) -o sbm sbm.o $(LDFLAGS)

clean:
	rm -f *.d
	rm -f *.o
	rm -f *~
	rm -f el2fg fg2fm fg2crs fg_lcc csr2fg sbm
