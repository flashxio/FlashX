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

include ../Makefile.common

CXXFLAGS += -I.. -I../include -I../libcommon
CFLAGS = -g -O3 -I.. -I../include -I../libcommon $(TRACE_FLAGS)
LDFLAGS := -L../libsafs -lsafs -L../libcommon -lcommon $(LDFLAGS)
LIBFILE = ../libsafs/libsafs.a ../libcommon/libcommon.a

TARGETS = SAFS-util

all: $(TARGETS)

SAFS-util: SAFS-util.o $(LIBFILE)
	$(CXX) -o SAFS-util SAFS-util.o $(LDFLAGS)

clean:
	rm -f *.o
	rm -f *~
	rm -f *.d
	rm -f $(TARGETS)
