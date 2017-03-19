# Copyright 2014 Open Connectome Project (http;//openconnecto.me)
# Written by Da Zheng (zhengda1936@gmail.com)
#
# This file is part of SAFSlib.
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

CFLAGS = -g -I../ $(TRACE_FLAGS)
LDFLAGS := -L.. -lsafs $(LDFLAGS)

UNITTEST = file_mapper_unit_test slab_allocator_test test_mem_tracker native_file_unit_test	\
		   safs_file_unit_test test_open_close test-io test-NUMA_buffer
CPPFLAGS := -MD
CXXFLAGS = -I.. -I../ -g -std=c++0x
SOURCE := $(wildcard *.c) $(wildcard *.cpp)
OBJS := $(patsubst %.c,%.o,$(patsubst %.cpp,%.o,$(SOURCE)))
DEPS := $(patsubst %.o,%.d,$(OBJS))
LIBFILE = ../libsafs.a

all: $(UNITTEST)

SAHT_unit_test: SAHT_unit_test.o $(LIBFILE)
	$(CXX) -o SAHT_unit_test SAHT_unit_test.o $(LDFLAGS)

GClock_unit_test: GClock_unit_test.o $(LIBFILE)
	$(CXX) -o GClock_unit_test GClock_unit_test.o $(LDFLAGS)

SA_expand_shrink_test: SA_expand_shrink_test.o $(LIBFILE)
	$(CXX) -o SA_expand_shrink_test SA_expand_shrink_test.o $(LDFLAGS)

slab_allocator_test: slab_allocator_test.o $(LIBFILE)
	$(CXX) -o slab_allocator_test slab_allocator_test.o $(LDFLAGS)

file_mapper_unit_test: file_mapper_unit_test.o $(LIBFILE)
	$(CXX) -o file_mapper_unit_test file_mapper_unit_test.o $(LDFLAGS)

test_mem_tracker: test_mem_tracker.o $(LIBFILE)
	$(CXX) -o test_mem_tracker test_mem_tracker.o $(LDFLAGS)

native_file_unit_test: native_file_unit_test.o $(LIBFILE)
	$(CXX) -o native_file_unit_test native_file_unit_test.o $(LDFLAGS)

safs_file_unit_test: safs_file_unit_test.o $(LIBFILE)
	$(CXX) -o safs_file_unit_test safs_file_unit_test.o $(LDFLAGS)

KV_store_unit_test: KV_store_unit_test.o $(LIBFILE)
	$(CXX) -o KV_store_unit_test KV_store_unit_test.o $(LDFLAGS)

test_open_close: test_open_close.o $(LIBFILE)
	$(CXX) -o test_open_close test_open_close.o $(LDFLAGS)

test-io: test-io.o $(LIBFILE)
	$(CXX) -o test-io test-io.o $(LDFLAGS)

test-NUMA_buffer: test-NUMA_buffer.o $(LIBFILE)
	$(CXX) -o test-NUMA_buffer test-NUMA_buffer.o $(LDFLAGS)

test:
	./slab_allocator_test
	./file_mapper_unit_test
	./test_mem_tracker
	./native_file_unit_test
	./test-NUMA_buffer
	mkdir -p /tmp/safs_data
	./safs_file_unit_test data_files.txt
	./test_open_close data_files.txt
	./test-io run_test.txt
	rm -R /tmp/safs_data

clean:
	rm -f *.o
	rm -f *.d
	rm -f *~
	rm -f $(UNITTEST)

-include $(DEPS) 
