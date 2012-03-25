CFLAGS = -g -O3 -DSTATISTICS #-DPROFILER
CXXFLAGS = -g -Wall -O3 -DSTATISTICS -DNUM_NODES=8 #-DPROFILER
LDFLAGS = -lpthread -lnuma #-laio #-lprofiler
HEADERS = aio_private.h direct_private.h global_cached_private.h mmap_private.h part_cached_private.h part_global_cached_private.h read_private.h thread_private.h workload.h rand_buf.h
CC = g++

all: rand-read rand-memcpy test create_file

create_file: create_file.o
	g++ -o create_file create_file.o -lnuma

rand-read.o: rand-read.cc cache.h tree_cache.h associative_cache.h cuckoo_cache.h wpaio.h LRU2Q.h $(HEADERS)
wpaio.o: wpaio.h

rand-read: rand-read.o wpaio.o
	g++ -o rand-read wpaio.o rand-read.o $(LDFLAGS)
rand-memcpy: rand-memcpy.c
	gcc -o rand-memcpy rand-memcpy.c $(LDFLAGS)

test: test.cc cuckoo_cache.h cache.h associative_cache.h tree_cache.h $(HEADERS)
	g++ -o test test.cc -g -lpthread

clean:
	rm *.o
