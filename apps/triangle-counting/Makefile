CFLAGS = -g -O3 -DSTATISTICS #-DPROFILER
ifdef MEMCHECK
TRACE_FLAGS = -faddress-sanitizer
else
	OMP_FLAG = -fopenmp
endif
TRACE_FLAGS += -fno-omit-frame-pointer # for better stack traces in error messages
TRACE_FLAGS += -fno-optimize-sibling-calls # disable tail call elimination
CLANG_FLAGS = -Wno-attributes
LDFLAGS = -lpthread $(TRACE_FLAGS) -lprofiler -rdynamic -L../../libcache -lcache -L.. -lgraph -laio -lnuma $(OMP_FLAG)
CXXFLAGS = -g -O3 -I../../include -I.. -I. -Wall -std=c++0x $(TRACE_FLAGS) $(CLANG_FLAGS) $(OMP_FLAG)
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

all: triangle-counting

triangle-counting: $(OBJS) ../libgraph.a
	$(CXX) -o triangle-counting $(OBJS) $(LDFLAGS)

build_lib:
	make -C libcache

unit_test: $(OBJS)
	make -C test

clean:
	rm -f *.d
	rm -f *.o
	rm -f *~
	rm -f triangle-counting

-include $(DEPS) 
