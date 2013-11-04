# Copyright 2013 Da Zheng
#
# This file is part of SAFSlib.
#
# SAFSlib is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# SAFSlib is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with SAFSlib.  If not, see <http://www.gnu.org/licenses/>.

ifdef MEMCHECK
TRACE_FLAGS = -faddress-sanitizer
CC = clang
CXX = clang++
else
CC = gcc
CXX = g++
endif
TRACE_FLAGS += -fno-omit-frame-pointer # for better stack traces in error messages
TRACE_FLAGS += -fno-optimize-sibling-calls # disable tail call elimination
CXXFLAGS = -g -O3 -I.. -I../include -std=c++0x $(TRACE_FLAGS)
CFLAGS = -g -O3 -I.. -I../include $(TRACE_FLAGS)

LDFLAGS = -L../libsafs -lsafs -lpthread -lnuma -laio $(TRACE_FLAGS) -lprofiler -lrt

TARGETS = SAFS-util

all: $(TARGETS)

SAFS-util: SAFS-util.o ../libsafs/libsafs.a
	$(CXX) -o SAFS-util SAFS-util.o $(LDFLAGS)

clean:
	rm -f *.o
	rm -f *~
	rm -f $(TARGETS)
