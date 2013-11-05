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
UNITTEST = file_mapper_unit_test slab_allocator_test test_mem_tracker native_file_unit_test	\
		   safs_file_unit_test unique_ptr_unit_test
CPPFLAGS := -MD
CXXFLAGS = -I.. -I../include -g -std=c++0x
SOURCE := $(wildcard *.c) $(wildcard *.cpp)
OBJS := $(patsubst %.c,%.o,$(patsubst %.cpp,%.o,$(SOURCE)))
DEPS := $(patsubst %.o,%.d,$(OBJS))
LIBFILE = ../libsafs/libsafs.a

all: $(UNITTEST)

SAHT_unit_test: SAHT_unit_test.o $(LIBFILE)
	$(CXX) -o SAHT_unit_test SAHT_unit_test.o $(LDFLAGS)

GClock_unit_test: GClock_unit_test.o $(LIBFILE)
	$(CXX) -o GClock_unit_test GClock_unit_test.o $(LDFLAGS)

SA_expand_shrink_test: SA_expand_shrink_test.o $(LIBFILE)
	$(CXX) -o SA_expand_shrink_test SA_expand_shrink_test.o $(LDFLAGS)

slab_allocator_test: slab_allocator_test.o $(LIBFILE)
	$(CXX) -o slab_allocator_test slab_allocator_test.o $(LDFLAGS)

file_mapper_unit_test: file_mapper_unit_test.o $(LIBFILE)
	$(CXX) -o file_mapper_unit_test file_mapper_unit_test.o $(LDFLAGS)

test_mem_tracker: test_mem_tracker.o $(LIBFILE)
	$(CXX) -o test_mem_tracker test_mem_tracker.o $(LDFLAGS)

native_file_unit_test: native_file_unit_test.o $(LIBFILE)
	$(CXX) -o native_file_unit_test native_file_unit_test.o $(LDFLAGS)

safs_file_unit_test: safs_file_unit_test.o $(LIBFILE)
	$(CXX) -o safs_file_unit_test safs_file_unit_test.o $(LDFLAGS)

unique_ptr_unit_test: unique_ptr_unit_test.o $(LIBFILE)
	$(CXX) -o unique_ptr_unit_test unique_ptr_unit_test.o $(LDFLAGS)

clean:
	rm -f *.o
	rm -f *.d
	rm -f *~
	rm -f $(UNITTEST)

-include $(DEPS) 
