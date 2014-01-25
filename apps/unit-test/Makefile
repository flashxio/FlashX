CC = gcc
CXX = g++

ifdef MEMCHECK
TRACE_FLAGS = -faddress-sanitizer
CC = clang
CXX = clang++
TRACE_FLAGS += -fno-omit-frame-pointer # for better stack traces in error messages
TRACE_FLAGS += -fno-optimize-sibling-calls # disable tail call elimination
CLANG_FLAGS = -Wno-attributes
endif

LDFLAGS = -L../../libsafs -lsafs -L.. -lgraph -lpthread -lnuma -laio -lprofiler $(TRACE_FLAGS) -lrt -rdynamic
UNITTEST = test-bitmap
CPPFLAGS := -MD
CXXFLAGS = -I.. -I../../include -g -std=c++0x
SOURCE := $(wildcard *.c) $(wildcard *.cpp)
OBJS := $(patsubst %.c,%.o,$(patsubst %.cpp,%.o,$(SOURCE)))
DEPS := $(patsubst %.o,%.d,$(OBJS))

all: $(UNITTEST)

test-bitmap: test-bitmap.o ../libgraph.a
	$(CXX) -o test-bitmap test-bitmap.o $(LDFLAGS)

clean:
	rm -f *.o
	rm -f *.d
	rm -f *~
	rm -f $(UNITTEST)

-include $(DEPS) 
