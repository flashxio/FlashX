CFLAGS = -g -O3
CXXFLAGS = -g -O3

all: rand-read rand-memcpy test create_file

create_file: create_file.o
	gcc -o create_file create_file.o

rand-read.o: rand-read.cc cache.h tree_cache.h associative_cache.h

rand-read: rand-read.o
	g++ -o rand-read rand-read.o -lpthread
rand-memcpy: rand-memcpy.o
	gcc -o rand-memcpy rand-memcpy.o -lpthread

test.o: test.cc cache.h

test: test.o
	g++ -o test test.o -lpthread
