include ../default.mk

###############################

DIR = ../

###############################

MKL_SEC = -I$(MKLROOT)/include -mkl=sequential 
MKL_PAR = -I$(MKLROOT)/include -mkl=parallel 

###############################

INCLUDE_SHARED = -I$(DIR)/Shared
LIBRARY_SHARED = -L$(DIR)/Shared -lclock -lsparsenew -lvector

###############################

INCLUDE_MPI = -I$(INC_MPI)
LIBRARY_MPI = -L$(LIB_MPI) -lmpi

################################

#CC      = mpicc -fopenmp
#F77     = mpifort -fopenmp
ifeq ($(CC),mpicc)
        OPENMP=-fopenmp
else
        OPENMP=-qopenmp
endif

CFLAGS  = $(OPENMP) -Wall -g -fPIC -fno-stack-protector -D_ILU0_ -D_PARDISO_MATCHING_ -D__UNDERSCORE__ -D_DOUBLE_REAL_ $(INCLUDE_SHARED) $(INC_MPI)
CLINKER = $(CC) 
FLINKER = $(F77) 
LDFLAGS = $(OPENMP)

LIBLIST = $(LIBRARY_SHARED) $(LIBRARY_MPI)

#######################################################################################################

EXECS =  Test 

default: Test 


SRC=  Lists.c SortFunction.c EliminationTree.c ILU0Factor.c ToolsILU0.c ToolsMPI.c ToolsMPI_OPENMP.c TreeOrderings.c SPDfactorMPI_OPENMP.c SPDsolverMPI_OPENMP.c Test_MPI_OPENMP.c

OBJ = $(SRC:.c=.o)

all: $(EXECS)

COMP:$(OBJ)

Test: Lists.o SortFunction.o EliminationTree.o ILU0Factor.o ToolsILU0.o ToolsMPI.o ToolsMPI_OPENMP.o TreeOrderings.o SPDfactorMPI_OPENMP.o SPDsolverMPI_OPENMP.o Test_MPI_OPENMP.o  
	$(FLINKER) $(LDFLAGS) $(mpicc -showme:link) $^ -o $@ $(LIBLIST)  

clean:
	/bin/rm -rf core *.o *~ PI* $(EXECS) aux[1,2,3]_* rank_* Trace[0-9]*

.c.o:
	echo compilando
	 $(CC) $(CFLAGS) $(mpicc -showme:compile) -c $<	

.f.o:
	echo compilando
	$(F77) $(FFLAGS) -c $<  -o $@

.F.o:
	echo compilando
	$(F77) $(FFLAGS) -c $<  -o $@

.h.h:
	echo compilando

