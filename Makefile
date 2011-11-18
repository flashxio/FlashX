CFLAGS = -g -O3 -DSTATISTICS #-DNDEBUG
CXXFLAGS = -g -O3 -DSTATISTICS

all: rand-read rand-memcpy test create_file

create_file: create_file.o
	gcc -o create_file create_file.o

rand-read.o: rand-read.cc cache.h tree_cache.h associative_cache.h cuckoo_cache.h

rand-read: rand-read.o
	g++ -o rand-read rand-read.o -lpthread
rand-memcpy: rand-memcpy.o
	gcc -o rand-memcpy rand-memcpy.o -lpthread

test: test.cc cuckoo_cache.h cache.h associative_cache.h tree_cache.h
	g++ -o test test.cc -g -lpthread

clean:
	rm *.o
