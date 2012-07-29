CC = gcc
CXX = g++

LDFLAGS = -lpthread -lnuma -laio

all: stat create_file rand-memcpy

stat: stat.o ../workload.h
	$(CXX) -o stat stat.o

create_file: create_file.o
	$(CXX) -o create_file create_file.o -lnuma

rand-memcpy: rand-memcpy.o
	$(CC) -o rand-memcpy rand-memcpy.o $(LDFLAGS)

clean:
	rm *.o
