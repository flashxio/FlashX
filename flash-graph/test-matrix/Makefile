# Copyright 2014 Open Connectome Project (http://openconnecto.me)
# Written by Disa Mhembere (disa@jhu.edu)
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
LDFLAGS := -L../matrix -lmatrix -L../libgraph-algs -L.. -lgraph -L../../libsafs -lsafs -L../../libcommon -lcommon -lrt $(OMP_FLAG) -lz $(LDFLAGS)
CXXFLAGS += -I../../include -I../../libcommon -I.. -I. $(OMP_FLAG) -I/usr/include/eigen3/ -DUSE_EIGEN

all: test-kmeans

test-kmeans: test-kmeans.o ../libgraph.a
	$(CXX) -o test-kmeans test-kmeans.o ../matrix/kmeans.o ../matrix/min-tri-kmeans.o\
		 ../matrix/kmeans_task_coordinator.o ../matrix/kmeans_coordinator.o \
		 ../matrix/base_kmeans_coordinator.o ../matrix/kmeans_task_thread.o \
		../matrix/kmeans_thread.o ../libgraph-algs/clusters.o \
		../libgraph-algs/dist_matrix.o $(LDFLAGS)

clean:
	rm -f *.d
	rm -f *.o
	rm -f *~
	rm -f test-kmeans

-include $(DEPS) 
