# Copyright 2016 neurodata (http://neurodata.io/)
# Written by Disa Mhembere (disa@jhu.edu)
#
# This file is part of k-par-means.
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

include ../../../Makefile.common

LDFLAGS := -L../../../../libkcommon -lcommon -L.. -lgraph \
	-L../../libsafs -lsafs -lrt -lboost_log $(LDFLAGS)
CXXFLAGS += -I../../../../libkcommon -I../../libsafs -I.. -I.

all: read_rrow_data

test: all
	./read_rrow_data 50 5 ../../../../../test-data/matrix_r50_c5_rrw.bin

read_rrow_data: read_rrow_data.o
	$(CXX) -o read_rrow_data read_rrow_data.o

clean:
	rm -f *.d
	rm -f *.o
	rm -f *~
	rm -f read_rrow_data

-include $(DEPS)
