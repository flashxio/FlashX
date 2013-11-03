# Copyright 2013 Da Zheng
#
# This file is part of SA-GraphLib.
#
# SA-GraphLib is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# SA-GraphLib is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with SA-GraphLib.  If not, see <http://www.gnu.org/licenses/>.

CFLAGS = -g -O3 -DSTATISTICS -DENABLE_AIO #-DPROFILER
ifdef MEMCHECK
TRACE_FLAGS = -faddress-sanitizer
else
	OMP_FLAG = -fopenmp
endif
TRACE_FLAGS += -fno-omit-frame-pointer # for better stack traces in error messages
TRACE_FLAGS += -fno-optimize-sibling-calls # disable tail call elimination
CLANG_FLAGS = -Wno-attributes
LDFLAGS = -lpthread -lnuma -laio $(TRACE_FLAGS) -lprofiler -rdynamic $(OMP_FLAG)
CXXFLAGS = -g -O3 -I../include -I. -Wall -std=c++0x $(TRACE_FLAGS) $(CLANG_FLAGS) $(OMP_FLAG)
CPPFLAGS := -MD

SOURCE := $(wildcard *.c) $(wildcard *.cpp)
OBJS := $(patsubst %.c,%.o,$(patsubst %.cpp,%.o,$(SOURCE)))
DEPS := $(patsubst %.o,%.d,$(OBJS))
MISSING_DEPS := $(filter-out $(wildcard $(DEPS)),$(DEPS))
MISSING_DEPS_SOURCES := $(wildcard $(patsubst %.d,%.c,$(MISSING_DEPS)) $(patsubst %.d,%.cc,$(MISSING_DEPS)))
ifdef MEMCHECK
CXXFLAGS += -DMEMCHECK
CC = clang
CXX = clang++
else
CC = gcc
CXX = g++
endif

all: libgraph graph-bfs triangle-counting tools

libgraph: $(OBJS)
	rm -f libgraph.a
	ar -cvq libgraph.a $(OBJS)

graph-bfs: libgraph
	$(MAKE) -C graph-bfs

triangle-counting: libgraph
	$(MAKE) -C triangle-counting

tools: libgraph
	$(MAKE) -C tools

clean:
	rm -f *.d
	rm -f *.o
	rm -f *~
	rm -f libgraph.a
	find -name core -delete
	$(MAKE) clean -C tools
	$(MAKE) clean -C graph-bfs
	$(MAKE) clean -C triangle-counting

-include $(DEPS) 
