include ../../Makefile.common

all: libeigen

CXXFLAGS += -I.. -I../../libsafs -I../../flash-graph

libeigen: $(OBJS)
	rm -f libeigen.a
	ar -cvq libeigen.a $(OBJS)

clean:
	rm -f *~
	rm -f *.o
	rm -f *.a
	rm -f *.d

-include $(DEPS) 
