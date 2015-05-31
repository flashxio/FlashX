include ../../Makefile.common

all: libeigen

OMP_FLAG = -fopenmp
CXXFLAGS += -I.. -I../../include -I../../libcommon -I../../flash-graph $(OMP_FLAG)

libeigen: $(OBJS)
	rm -f libeigen.a
	ar -cvq libeigen.a $(OBJS)

clean:
	rm *~
	rm *.o
	rm *.a

-include $(DEPS) 
