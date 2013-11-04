ifdef MEMCHECK
TRACE_FLAGS = -faddress-sanitizer
CC = clang
CXX = clang++
else
CC = gcc
CXX = g++
endif
TRACE_FLAGS += -fno-omit-frame-pointer # for better stack traces in error messages
TRACE_FLAGS += -fno-optimize-sibling-calls # disable tail call elimination
CXXFLAGS = -g -O3 -I.. -I../include -std=c++0x $(TRACE_FLAGS)
CFLAGS = -g -O3 -I.. -I../include $(TRACE_FLAGS)

LDFLAGS = -L../libcache -lcache -lpthread -lnuma -laio $(TRACE_FLAGS) -lprofiler -lrt

TARGETS = SAFS-util

all: $(TARGETS)

SAFS-util: SAFS-util.o ../libcache/libcache.a
	$(CXX) -o SAFS-util SAFS-util.o $(LDFLAGS)

clean:
	rm -f *.o
	rm -f *~
	rm -f $(TARGETS)
