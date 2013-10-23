CFLAGS = -g -O3 -DSTATISTICS #-DPROFILER
ifdef MEMCHECK
TRACE_FLAGS = -faddress-sanitizer
endif
TRACE_FLAGS += -fno-omit-frame-pointer # for better stack traces in error messages
TRACE_FLAGS += -fno-optimize-sibling-calls # disable tail call elimination
CLANG_FLAGS = -Wno-attributes
LDFLAGS = -lpthread $(TRACE_FLAGS) -lprofiler -rdynamic -L../../libcache -lcache -L.. -lgraph -laio -lnuma
CXXFLAGS = -g -O3 -I../../include -I.. -I. -Wall -std=c++0x $(TRACE_FLAGS) $(CLANG_FLAGS) -DNUM_NODES=1 -DNCPUS=0 -DPROFILER -DSTATISTICS
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

all: graph-bfs

graph-bfs: $(OBJS) ../libgraph.a
	$(CXX) -o graph-bfs $(OBJS) $(LDFLAGS)

build_lib:
	make -C libcache

unit_test: $(OBJS)
	make -C test

clean:
	rm -f *.d
	rm -f *.o
	rm -f *~
	rm graph-bfs

-include $(DEPS) 
