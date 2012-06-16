CFLAGS = -g -O3 -DSTATISTICS -DENABLE_AIO #-DPROFILER
#TRACE_FLAGS = -faddress-sanitizer 
#TRACE_FLAGS += -fno-omit-frame-pointer # for better stack traces in error messages
#TRACE_FLAGS += -fno-optimize-sibling-calls # disable tail call elimination
CLANG_FLAGS = -Wno-attributes
LDFLAGS = -lpthread -lnuma -laio $(TRACE_FLAGS) #-lprofiler
CXXFLAGS = -g -Inbds.0.4.3/include/ -Wall $(TRACE_FLAGS) $(CLANG_FLAGS) -DSTATISTICS -DENABLE_AIO -DNUM_NODES=1 -DNCPUS=0 #-DPROFILER
HEADERS = aio_private.h direct_private.h global_cached_private.h part_global_cached_private.h read_private.h thread_private.h workload.h rand_buf.h
CC = clang
CXX = clang++
OBJECTS = part_global_cached_private.o cuckoo_cache.o associative_cache.o workload.o global_cached_private.o direct_private.o read_private.o rand_buf.o memory_manager.o wpaio.o rand-read.o aio_private.o hash_index_cache.o messaging.o disk_read_thread.o thread_private.o

all: rand-read rand-memcpy create_file

create_file: create_file.o
	$(CXX) -o create_file create_file.o -lnuma

rand-read.o: rand-read.cc cache.h tree_cache.h associative_cache.h cuckoo_cache.h wpaio.h LRU2Q.h $(HEADERS)
wpaio.o: wpaio.h

rand-read: $(OBJECTS)
	$(CXX) -o rand-read -static $(OBJECTS) $(LDFLAGS) -Lnbds.0.4.3/ -lnbds
rand-memcpy: rand-memcpy.c
	$(CC) -o rand-memcpy rand-memcpy.c $(LDFLAGS)

test: test.cc cuckoo_cache.h cache.h associative_cache.h tree_cache.h $(HEADERS)
	$(CXX) -o test test.cc -g -lpthread -DNUM_NODES=1

clean:
	rm *.o
