CFLAGS = -g -O3 -DSTATISTICS #-DPROFILER
ifdef MEMCHECK
TRACE_FLAGS = -faddress-sanitizer
endif
TRACE_FLAGS += -fno-omit-frame-pointer # for better stack traces in error messages
TRACE_FLAGS += -fno-optimize-sibling-calls # disable tail call elimination
CLANG_FLAGS = -Wno-attributes
LDFLAGS = -lpthread $(TRACE_FLAGS) -lprofiler -rdynamic -Llibcache -lcache -laio -lnuma -lrt
CXXFLAGS = -g -O3 -Iinclude -I. -Wall -std=c++0x $(TRACE_FLAGS) $(CLANG_FLAGS) -DPROFILER -DSTATISTICS
ifdef USE_NBDS
LDFLAGS += -Lnbds.0.4.3/ -lnbds
CXXFLAGS += -Inbds.0.4.3/include/ -DUSE_NBDS
endif
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

all: rand-read unit_test tools apps

rand-read: $(OBJS) build_lib
	$(CXX) -o rand-read $(OBJS) $(LDFLAGS)

build_lib:
	$(MAKE) -C libcache

unit_test: build_lib
ifndef MEMCHECK
	$(MAKE) -C test
endif

tools: build_lib
	$(MAKE) -C tools

apps: build_lib
	$(MAKE) -C apps

clean:
	rm -f *.d
	rm -f *.o
	rm -f *~
	rm -f include/*~
	find -name core -delete
	make --ignore-errors -C test clean
	make --ignore-errors -C libcache clean
	make --ignore-errors -C tools clean
	make --ignore-errors -C apps clean

-include $(DEPS) 
