include ../../Makefile.common

OMP_FLAG = -fopenmp
LDFLAGS := -L.. -lgraph -L../../libsafs -lsafs -lrt $(OMP_FLAG) $(LDFLAGS)
CXXFLAGS = -I.. -I../../include -g -std=c++0x

SOURCE := $(wildcard *.c) $(wildcard *.cpp)
OBJS := $(patsubst %.c,%.o,$(patsubst %.cpp,%.o,$(SOURCE)))
DEPS := $(patsubst %.o,%.d,$(OBJS))

UNITTEST = test-bitmap test-partitioner

all: $(UNITTEST)

test-bitmap: test-bitmap.o ../libgraph.a
	$(CXX) -o test-bitmap test-bitmap.o $(LDFLAGS)

test-partitioner: test-partitioner.o ../libgraph.a
	$(CXX) -o test-partitioner test-partitioner.o $(LDFLAGS)

clean:
	rm -f *.o
	rm -f *.d
	rm -f *~
	rm -f $(UNITTEST)

-include $(DEPS) 
