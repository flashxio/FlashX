CFLAGS = -g -O3 -DSTATISTICS -DENABLE_AIO #-DPROFILER
ifdef MEMCHECK
TRACE_FLAGS = -faddress-sanitizer
endif
TRACE_FLAGS += -fno-omit-frame-pointer # for better stack traces in error messages
TRACE_FLAGS += -fno-optimize-sibling-calls # disable tail call elimination
CLANG_FLAGS = -Wno-attributes
LDFLAGS = -lpthread -lnuma -laio $(TRACE_FLAGS) -lprofiler -rdynamic
CXXFLAGS = -g -O3 -I../include -I. -I../nbds.0.4.3/include/ -Wall -std=c++0x $(TRACE_FLAGS) $(CLANG_FLAGS) -DSTATISTICS -DENABLE_AIO -DNUM_NODES=1 -DNCPUS=0 -DPROFILER
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

all: libcache

libcache: $(OBJS)
	ar -cvq libcache.a $(OBJS)

clean:
	rm -f *.d
	rm -f *.o
	rm -f *~
	rm libcache.a

-include $(DEPS) 
