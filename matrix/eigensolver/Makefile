include ../../Makefile.common

all: BlockDavidsonEpetra

OMP_FLAG = -fopenmp
CXXFLAGS += -I.. -I../../include -I../../libcommon -I../../flash-graph $(OMP_FLAG)

LDFLAGS := -lteuchoscomm -lteuchoscore -lteuchosparameterlist -lteuchoskokkoscomm	\
		 -lteuchosnumerics -lanasazi -lanasaziepetra -lblas -llapack -lkokkoscore	\
		 $(OMP_FLAG) -L.. -lFMatrix -L../../flash-graph -lgraph -L../../libsafs -lsafs	\
		 -L../../libcommon -lcommon $(LDFLAGS)
LDFLAGS += -lstxxl -lz -lnuma -laio

BlockDavidsonEpetra: BlockDavidsonEpetra.o
	$(CC) BlockDavidsonEpetra.o -o BlockDavidsonEpetra $(LDFLAGS)

clean:
	rm *~
	rm *.o
	rm BlockDavidsonEpetra

-include $(DEPS) 
