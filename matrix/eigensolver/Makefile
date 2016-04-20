include ../../Makefile.common

all: libeigen

OMP_FLAG = -fopenmp
CXXFLAGS += -I.. -I../../libsafs -I../../flash-graph $(OMP_FLAG)

libeigen: $(OBJS)
	rm -f libeigen.a
	ar -cvq libeigen.a $(OBJS)

clean:
	rm -f *~
	rm -f *.o
	rm -f *.a
	rm -f *.d

-include $(DEPS) 
