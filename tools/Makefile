CFLAGS = -g
CC = gcc
CXX = g++

LDFLAGS = -lpthread -lnuma -laio

all: stat create_file rand-memcpy rand-memcpy-agg

stat: stat.o ../workload.h
	$(CXX) -o stat stat.o

create_file: create_file.o
	$(CXX) -o create_file create_file.o -lnuma

rand-memcpy: rand-memcpy.o
	$(CC) -o rand-memcpy rand-memcpy.o $(LDFLAGS)

rand-memcpy-agg: rand-memcpy-agg.o
	$(CC) -o rand-memcpy-agg rand-memcpy-agg.o $(LDFLAGS)

clean:
	rm *.o
