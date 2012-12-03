CFLAGS = -g -O3 -DSTATISTICS -DENABLE_AIO #-DPROFILER
ifdef MEMCHECK
TRACE_FLAGS = -faddress-sanitizer
endif
TRACE_FLAGS += -fno-omit-frame-pointer # for better stack traces in error messages
TRACE_FLAGS += -fno-optimize-sibling-calls # disable tail call elimination
CLANG_FLAGS = -Wno-attributes
LDFLAGS = -lpthread -lnuma -laio $(TRACE_FLAGS) -lprofiler
CXXFLAGS = -g -Inbds.0.4.3/include/ -Wall $(TRACE_FLAGS) $(CLANG_FLAGS) -DSTATISTICS -DENABLE_AIO -DNUM_NODES=1 -DNCPUS=0 -DPROFILER
CPPFLAGS := -MD

SOURCE := $(wildcard *.c) $(wildcard *.cpp)
OBJS := $(patsubst %.c,%.o,$(patsubst %.cpp,%.o,$(SOURCE)))
DEPS := $(patsubst %.o,%.d,$(OBJS))
MISSING_DEPS := $(filter-out $(wildcard $(DEPS)),$(DEPS))
MISSING_DEPS_SOURCES := $(wildcard $(patsubst %.d,%.c,$(MISSING_DEPS)) $(patsubst %.d,%.cc,$(MISSING_DEPS)))
ifdef MEMCHECK
CC = clang
CXX = clang++
else
CC = gcc
CXX = g++
endif

all: rand-read

rand-read: $(OBJS)
	$(CXX) -o rand-read $(OBJS) $(LDFLAGS) -Lnbds.0.4.3/ -lnbds

unit_test: $(OBJS)
	make -C test

clean:
	rm -f *.d
	rm -f *.o
	rm -f *~
	make -C test clean

-include $(DEPS) 
