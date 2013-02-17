CXXFLAGS = -g -std=c++0x -O2 -I../..
CC = gcc
CXX = g++

all: BeladyAlgo

BeladyAlgo: BeladyAlgo.o BeladyAlgoFile.o BeladyAlgoTest.o
	$(CXX) -o BeladyAlgo BeladyAlgo.o BeladyAlgoFile.o -std=c++0x -lprofiler
	$(CXX) -o BeladyAlgoTest BeladyAlgo.o BeladyAlgoTest.o -std=c++0x

clean:
	rm *.o
