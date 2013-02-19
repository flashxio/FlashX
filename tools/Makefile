CFLAGS = -g
CXXFLAGS = -g -I.. -I../nbds.0.4.3/include/
CC = gcc
CXX = g++

LDFLAGS = -lpthread -lnuma -laio

all: stat create_file rand-memcpy rand-memcpy-agg workload-gen memory-fill cache_evaluator

stat: stat.o ../workload.h
	$(CXX) -o stat stat.o

create_file: create_file.o
	$(CXX) -o create_file create_file.o -lnuma

rand-memcpy: rand-memcpy.o
	$(CC) -o rand-memcpy rand-memcpy.o $(LDFLAGS)

rand-memcpy-agg: rand-memcpy-agg.o
	$(CC) -o rand-memcpy-agg rand-memcpy-agg.o $(LDFLAGS)

workload-gen: workload-gen.o
	$(CXX) -o workload-gen workload-gen.o

memory-fill: memory-fill.o
	$(CXX) -o memory-fill memory-fill.o

cache_evaluator: cache_evaluator.o ../associative_cache.o ../workload.o ../slab_allocator.o ../hash_index_cache.o ../memory_manager.o ../messaging.o ../gclock.o
	$(CXX) -o cache_evaluator cache_evaluator.o ../associative_cache.o ../workload.o ../slab_allocator.o ../hash_index_cache.o ../memory_manager.o ../messaging.o ../gclock.o $(LDFLAGS) -L../nbds.0.4.3/ -lnbds

clean:
	rm *.o
