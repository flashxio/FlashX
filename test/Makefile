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

CFLAGS = -g -I../include $(TRACE_FLAGS)

LDFLAGS = -L../libsafs -lsafs -lpthread -lnuma -laio -lprofiler $(TRACE_FLAGS) -lrt -rdynamic
CPPFLAGS := -MD
CXXFLAGS = -I.. -I../include -g
SOURCE := $(wildcard *.c) $(wildcard *.cpp)
OBJS := $(patsubst %.c,%.o,$(patsubst %.cpp,%.o,$(SOURCE)))
DEPS := $(patsubst %.o,%.d,$(OBJS))

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
