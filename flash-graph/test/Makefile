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

LDFLAGS := -L../../../../libcommon -lcommon -L.. -lgraph \
	-L../../libsafs -lsafs -lrt -lz -lboost_log $(LDFLAGS)
CXXFLAGS += -I../../../../libcommon -I../../libsafs -I.. -I.

TESTFILES = test_load_balancer test_comm test_knors_index make_page_row_header\
			test_sem_kmeans test_partition_cache

all: $(TESTFILES)

test_load_balancer: test_load_balancer.o ../libgraph.a
	$(CXX) -o test_load_balancer test_load_balancer.o $(LDFLAGS)

test_comm: test_comm.o ../libgraph.a
	$(CXX) -o test_comm test_comm.o $(LDFLAGS)

test_knors_index: test_knors_index.o
	$(CXX) -o test_knors_index test_knors_index.o $(LDFLAGS)

test_partition_cache: test_partition_cache.o
	$(CXX) -o test_partition_cache test_partition_cache.o $(LDFLAGS)

make_page_row_header: make_page_row_header.o
	$(CXX) -o make_page_row_header make_page_row_header.o $(LDFLAGS)

test_sem_kmeans: test_sem_kmeans.o ../libgraph.a
	$(CXX) -o test_sem_kmeans test_sem_kmeans.o $(LDFLAGS)

clean:
	rm -f *.d
	rm -f *.o
	rm -f *~
	rm -f $(TESTFILES)

-include $(DEPS)
