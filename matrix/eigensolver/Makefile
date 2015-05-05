include ../../Makefile.common

all: BlockDavidsonTpetra BlockDavidsonFM

OMP_FLAG = -fopenmp
CXXFLAGS += -I.. -I../../include -I../../libcommon -I../../flash-graph $(OMP_FLAG)

LDFLAGS := -lepetra -ltpetra -lteuchoscomm -lteuchoskokkoscomm -lteuchosnumerics 	\
		 -lteuchosparameterlist -lteuchoscore -lanasazi -lanasaziepetra -lblas -llapack -lkokkoscore	\
		 -lkokkoscontainers -lteuchoskokkoscompat	\
		 $(OMP_FLAG) -L.. -lFMatrix -L../../flash-graph -lgraph -L../../libsafs -lsafs	\
		 -L../../libcommon -lcommon $(LDFLAGS)
LDFLAGS += -lstxxl -lz -lnuma -laio

BlockDavidsonTpetra: BlockDavidsonTpetra.o
	$(CXX) BlockDavidsonTpetra.o -o BlockDavidsonTpetra $(LDFLAGS)

BlockDavidsonFM: BlockDavidsonFM.o block_dense_matrix.o
	$(CXX) BlockDavidsonFM.o block_dense_matrix.o -o BlockDavidsonFM $(LDFLAGS)

clean:
	rm *~
	rm *.o
	rm BlockDavidsonTpetra
	rm BlockDavidsonFM

-include $(DEPS) 
