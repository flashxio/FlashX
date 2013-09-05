CC = gcc

ifdef MEMCHECK
TRACE_FLAGS = -faddress-sanitizer
CC = clang
CXX = clang++
TRACE_FLAGS += -fno-omit-frame-pointer # for better stack traces in error messages
TRACE_FLAGS += -fno-optimize-sibling-calls # disable tail call elimination
CLANG_FLAGS = -Wno-attributes
endif

CFLAGS = -g -I../include $(TRACE_FLAGS)

LDFLAGS = -L../libcache -lcache -lpthread -lnuma -laio -lprofiler -L../nbds.0.4.3/ -lnbds $(TRACE_FLAGS)
UNITTEST = SAHT_unit_test GClock_unit_test SA_expand_shrink_test file_mapper_unit_test slab_allocator_test test
CPPFLAGS := -MD
CXXFLAGS = -I.. -g
SOURCE := $(wildcard *.c) $(wildcard *.cpp)
OBJS := $(patsubst %.c,%.o,$(patsubst %.cpp,%.o,$(SOURCE)))
DEPS := $(patsubst %.o,%.d,$(OBJS))

test_c_interface: test_c_interface.o ../libcache/libcache.a
	$(CXX) -o test_c_interface test_c_interface.o $(LDFLAGS)

all: $(UNITTEST)

SAHT_unit_test: SAHT_unit_test.o ../SA_hash_table.o
	$(CXX) -o SAHT_unit_test SAHT_unit_test.o ../SA_hash_table.o $(LDFLAGS)

GClock_unit_test: GClock_unit_test.o ../gclock.o ../slab_allocator.o
	$(CXX) -o GClock_unit_test GClock_unit_test.o ../gclock.o ../slab_allocator.o $(LDFLAGS)

SA_expand_shrink_test: SA_expand_shrink_test.o ../associative_cache.o ../memory_manager.o ../slab_allocator.o ../messaging.o ../thread.o ../common.o
	$(CXX) -o SA_expand_shrink_test SA_expand_shrink_test.o ../associative_cache.o ../memory_manager.o ../slab_allocator.o ../messaging.o ../thread.o ../common.o $(LDFLAGS)

file_mapper_unit_test: file_mapper_unit_test.o
	$(CXX) -o file_mapper_unit_test file_mapper_unit_test.o

slab_allocator_test: slab_allocator_test.o
	$(CXX) -o slab_allocator_test slab_allocator_test.o $(LDFLAGS)

clean:
	rm *.o

-include $(DEPS) 
