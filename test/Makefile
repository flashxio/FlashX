# Copyright 2014 Open Connectome Project (http://openconnecto.me)
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

include ../Makefile.common

CFLAGS = -g -I../include $(TRACE_FLAGS)

LDFLAGS := -L../libsafs -lsafs $(LDFLAGS)
CXXFLAGS += -I.. -I../include

all: test_rand_io

test_rand_io: test_rand_io.o thread_private.o ../libsafs/libsafs.a
	$(CXX) -o test_rand_io test_rand_io.o thread_private.o $(LDFLAGS)

clean:
	rm -f *.o
	rm -f *.d
	rm -f *~
	rm -f test_rand_io

-include $(DEPS) 
