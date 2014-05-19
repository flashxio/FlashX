# Copyright 2014 Open Connectome Project (http://openconnecto.me)
# Written by Disa Mhembere (disa@jhu.edu)
#            Da Zheng (zhengda1936@gmail.com)
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
