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

CC      = mcc --ompss-v2 -keep --debug --no-copy-deps 
F77     = mfc --ompss-v2 -keep --debug --no-copy-deps 

CFLAGS  = -Wall -g -fPIC -fno-stack-protector -D_ILU0_ -D_PARDISO_MATCHING_ -D__UNDERSCORE__ -D_DOUBLE_REAL_ $(INCLUDE_SHARED) -I$(INC_MPI) 
CLINKER = mcc --ompss-v2 -keep --debug 
FLINKER = mfc --ompss-v2 -keep --debug 
LDFLAGS = --Wl,-nofor-main --no-copy-deps
LDFLAGS = --no-copy-deps

LIBLIST = $(LIBRARY_SHARED) $(LIBRARY_MPI)

#######################################################################################################

EXECS =  Test 

default: Test 


SRC=  Lists.c EliminationTree.c ILU0Factor.c ToolsILU0.c ToolsMPI.c ToolsMPI_OMPSS.c TreeOrderings.c SPDfactorMPI_OMPSS.c SPDsolverMPI_OMPSS.c Test_MPI_OMPSS.c

OBJ = $(SRC:.c=.o)

all: $(EXECS)

COMP:$(OBJ)

Test: Lists.o EliminationTree.o ILU0Factor.o ToolsILU0.o ToolsMPI.o ToolsMPI_OMPSS.o TreeOrderings.o SPDfactorMPI_OMPSS.o SPDsolverMPI_OMPSS.o Test_MPI_OMPSS.o  
	$(FLINKER) $(LDFLAGS) $(OPTFLAGS) $(mpicc -showme:link) $^ -o $@ $(LIBLIST)  

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

