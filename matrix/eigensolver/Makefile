include ../../Makefile.common

all: BlockDavidsonEpetra BlockDavidsonFM

OMP_FLAG = -fopenmp
CXXFLAGS += -I.. -I../../include -I../../libcommon -I../../flash-graph $(OMP_FLAG)

LDFLAGS := -lepetra -ltpetra -lteuchoscomm -lteuchoskokkoscomm -lteuchosnumerics 	\
		 -lteuchosparameterlist -lteuchoscore -lanasazi -lanasaziepetra -lblas -llapack -lkokkoscore	\
		 -lkokkoscontainers -lteuchoskokkoscompat	\
		 $(OMP_FLAG) -L.. -lFMatrix -L../../flash-graph -lgraph -L../../libsafs -lsafs	\
		 -L../../libcommon -lcommon $(LDFLAGS)
LDFLAGS += -lstxxl -lz -lnuma -laio

BlockDavidsonEpetra: BlockDavidsonEpetra.o
	$(CXX) BlockDavidsonEpetra.o -o BlockDavidsonEpetra $(LDFLAGS)

BlockDavidsonFM: BlockDavidsonFM.o
	$(CXX) BlockDavidsonFM.o -o BlockDavidsonFM $(LDFLAGS)

clean:
	rm *~
	rm *.o
	rm BlockDavidsonEpetra
	rm BlockDavidsonFM

-include $(DEPS) 
