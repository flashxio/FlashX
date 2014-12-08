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
LDFLAGS := -L.. -lgraph -L../../libsafs -lsafs -L../../libcommon -lcommon -lrt $(OMP_FLAG) -lstxxl -lz $(LDFLAGS)
CXXFLAGS += -I../../include -I../../libcommon -I.. -I. $(OMP_FLAG)

all: test_load_balancer test_msg_passing test_get test_push

test_load_balancer: test_load_balancer.o ../libgraph.a
	$(CXX) -o test_load_balancer test_load_balancer.o $(LDFLAGS)

test_msg_passing: test_msg_passing.o ../libgraph.a
	$(CXX) -o test_msg_passing test_msg_passing.o $(LDFLAGS)

test_get: test_get.o ../libgraph.a
	$(CXX) -o test_get test_get.o $(LDFLAGS)

test_push: test_push.o ../libgraph.a
	$(CXX) -o test_push test_push.o $(LDFLAGS)

clean:
	rm -f *.d
	rm -f *.o
	rm -f *~
	rm -f test_load_balancer
	rm -f test_msg_passing
	rm -f test_get
	rm -f test_push

-include $(DEPS) 
