CFLAGS = -g -O3 -DSTATISTICS -DENABLE_AIO #-DPROFILER
CXXFLAGS = -g -O3 -Wall -DSTATISTICS -DENABLE_AIO -DNUM_NODES=1 -DNCPUS=12 -DCPU_AFFINITY #-DPROFILER
LDFLAGS = -lpthread -lnuma -laio #-lprofiler
HEADERS = aio_private.h direct_private.h global_cached_private.h mmap_private.h part_cached_private.h part_global_cached_private.h read_private.h thread_private.h workload.h rand_buf.h
CC = g++
OBJECTS = part_global_cached_private.o cuckoo_cache.o associative_cache.o workload.o global_cached_private.o direct_private.o read_private.o rand_buf.o memory_manager.o wpaio.o rand-read.o aio_private.o

all: rand-read rand-memcpy create_file

create_file: create_file.o
	g++ -o create_file create_file.o -lnuma

rand-read.o: rand-read.cc cache.h tree_cache.h associative_cache.h cuckoo_cache.h wpaio.h LRU2Q.h $(HEADERS)
wpaio.o: wpaio.h

rand-read: $(OBJECTS)
	g++ -o rand-read $(OBJECTS) $(LDFLAGS)
rand-memcpy: rand-memcpy.c
	gcc -o rand-memcpy rand-memcpy.c $(LDFLAGS)

test: test.cc cuckoo_cache.h cache.h associative_cache.h tree_cache.h $(HEADERS)
	g++ -o test test.cc -g -lpthread -DNUM_NODES=1

clean:
	rm *.o
