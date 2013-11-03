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

LDFLAGS = -L../libcache -lcache -lpthread -lnuma -laio -lprofiler -L../nbds.0.4.3/ -lnbds $(TRACE_FLAGS) -lrt -rdynamic
UNITTEST = file_mapper_unit_test slab_allocator_test test_mem_tracker
CPPFLAGS := -MD
CXXFLAGS = -I.. -I../include -g
SOURCE := $(wildcard *.c) $(wildcard *.cpp)
OBJS := $(patsubst %.c,%.o,$(patsubst %.cpp,%.o,$(SOURCE)))
DEPS := $(patsubst %.o,%.d,$(OBJS))

all: $(UNITTEST)

SAHT_unit_test: SAHT_unit_test.o ../SA_hash_table.o
	$(CXX) -o SAHT_unit_test SAHT_unit_test.o ../SA_hash_table.o $(LDFLAGS)

GClock_unit_test: GClock_unit_test.o
	$(CXX) -o GClock_unit_test GClock_unit_test.o $(LDFLAGS)

SA_expand_shrink_test: SA_expand_shrink_test.o
	$(CXX) -o SA_expand_shrink_test SA_expand_shrink_test.o $(LDFLAGS)

slab_allocator_test: slab_allocator_test.o
	$(CXX) -o slab_allocator_test slab_allocator_test.o $(LDFLAGS)

file_mapper_unit_test: file_mapper_unit_test.o
	$(CXX) -o file_mapper_unit_test file_mapper_unit_test.o $(LDFLAGS)

test_mem_tracker: test_mem_tracker.o ../libcache/libcache.a
	$(CXX) -o test_mem_tracker test_mem_tracker.o $(LDFLAGS)

clean:
	rm -f *.o
	rm -f *.d
	rm -f *~
	rm -f $(UNITTEST)

-include $(DEPS) 
