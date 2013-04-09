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
CXXFLAGS = -g -O3 -I.. -I../nbds.0.4.3/include/ -std=c++0x $(TRACE_FLAGS)
CFLAGS = -g -O3 -I.. $(TRACE_FLAGS)

LDFLAGS = -lpthread -lnuma -laio $(TRACE_FLAGS) -lprofiler

all: stat create_file create_files mem-read rand-memcpy rand-memcpy-agg workload-gen memory-fill cache_evaluator java_dump2c_dump eval_expand_SA_cache

stat: stat.o ../workload.h
	$(CXX) -o stat stat.o $(LDFLAGS)

create_file: create_file.o
	$(CXX) -o create_file create_file.o $(LDFLAGS)

create_files: create_files.o ../common.o ../thread.o
	$(CXX) -o create_files create_files.o ../common.o ../thread.o $(LDFLAGS)

mem-read: mem-read.o
	$(CC) -o mem-read mem-read.o $(LDFLAGS)

rand-memcpy: rand-memcpy.o
	$(CC) -o rand-memcpy rand-memcpy.o $(LDFLAGS)

rand-memcpy-agg: rand-memcpy-agg.o
	$(CC) -o rand-memcpy-agg rand-memcpy-agg.o $(LDFLAGS)

workload-gen: workload-gen.o
	$(CXX) -o workload-gen workload-gen.o $(LDFLAGS)

memory-fill: memory-fill.o
	$(CXX) -o memory-fill memory-fill.o $(LDFLAGS)

java_dump2c_dump: java_dump2c_dump.o ../workload.o ../common.o
	$(CXX) -o java_dump2c_dump java_dump2c_dump.o ../workload.o ../common.o $(LDFLAGS)

cache_evaluator: cache_evaluator.o ../associative_cache.o ../workload.o ../slab_allocator.o ../hash_index_cache.o ../memory_manager.o ../messaging.o ../gclock.o ../thread.o ../common.o
	$(CXX) -o cache_evaluator cache_evaluator.o ../associative_cache.o ../workload.o ../slab_allocator.o ../hash_index_cache.o ../memory_manager.o ../messaging.o ../gclock.o ../thread.o ../common.o $(LDFLAGS) -L../nbds.0.4.3/ -lnbds

eval_expand_SA_cache: eval_expand_SA_cache.o ../associative_cache.o ../workload.o ../slab_allocator.o ../memory_manager.o ../messaging.o ../thread.o ../common.o
	$(CXX) -o eval_expand_SA_cache eval_expand_SA_cache.o ../associative_cache.o ../workload.o ../slab_allocator.o ../memory_manager.o ../messaging.o ../thread.o ../common.o $(LDFLAGS)

clean:
	rm *.o
