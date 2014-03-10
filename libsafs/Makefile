#
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

CXXFLAGS += -I../include

all: libsafs

libsafs: $(OBJS)
	rm -f libsafs.a
	ar -cvq libsafs.a $(OBJS)

clean:
	rm -f *.d
	rm -f *.o
	rm -f *~
	rm libsafs.a

-include $(DEPS) 
