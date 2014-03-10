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

include ../Makefile.common

CFLAGS = -g -I../include $(TRACE_FLAGS)

LDFLAGS := -L../libsafs -lsafs $(LDFLAGS)
CXXFLAGS += -I.. -I../include

all: test_c_interface test_rand_io

test_c_interface: test_c_interface.o ../libsafs/libsafs.a
	$(CXX) -o test_c_interface test_c_interface.o $(LDFLAGS)

test_rand_io: test_rand_io.o thread_private.o ../libsafs/libsafs.a
	$(CXX) -o test_rand_io test_rand_io.o thread_private.o $(LDFLAGS)

clean:
	rm -f *.o
	rm -f *.d
	rm -f *~
	rm -f test_c_interface
	rm -f test_rand_io

-include $(DEPS) 
