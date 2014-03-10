include ../../Makefile.common

OMP_FLAG = -fopenmp
LDFLAGS := -L../../libsafs -lsafs -L.. -lgraph -lrt $(OMP_FLAG) $(LDFLAGS)
CXXFLAGS = -I.. -I../../include -g -std=c++0x

UNITTEST = test-bitmap
SOURCE := $(wildcard *.c) $(wildcard *.cpp)
OBJS := $(patsubst %.c,%.o,$(patsubst %.cpp,%.o,$(SOURCE)))
DEPS := $(patsubst %.o,%.d,$(OBJS))

all: $(UNITTEST)

test-bitmap: test-bitmap.o ../libgraph.a
	$(CXX) -o test-bitmap test-bitmap.o $(LDFLAGS)

clean:
	rm -f *.o
	rm -f *.d
	rm -f *~
	rm -f $(UNITTEST)

-include $(DEPS) 
