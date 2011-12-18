CFLAGS = -g -O3 -DSTATISTICS #-DPROFILER
CXXFLAGS = -g -O3 -DSTATISTICS #-DPROFILER
LDFLAGS = -lpthread -laio #-lprofiler
CC = g++

all: rand-read rand-memcpy test create_file

create_file: create_file.o
	gcc -o create_file create_file.o

rand-read.o: rand-read.cc cache.h tree_cache.h associative_cache.h cuckoo_cache.h wpaio.h
wpaio.o: wpaio.h

rand-read: rand-read.o wpaio.o
	g++ -o rand-read wpaio.o rand-read.o $(LDFLAGS)
rand-memcpy: rand-memcpy.c
	gcc -o rand-memcpy rand-memcpy.c $(LDFLAGS)

test: test.cc cuckoo_cache.h cache.h associative_cache.h tree_cache.h
	g++ -o test test.cc -g -lpthread

clean:
	rm *.o
