# Compilers
CXX=g++
MCXX=mcxx
MPICXX=mpicxx

# Use the Intel MPI compiler if possible
ifdef I_MPI_ROOT
MPICXX=mpiicpc
endif

# Set the default block size (elements)
BSX?=1024
BSY?=$(BSX)

# Preprocessor flags
CPPFLAGS=-Isrc -DBSX=$(BSX) -DBSY=$(BSY)

# Compiler flags
CFLAGS=-O3 -std=c++11
MCCFLAGS=--ompss-2 $(CFLAGS) --Wn,-O3,-std=c++11

# Linker flags
LDFLAGS=-lrt -lm

# Link to the multithreaded MPI
MPI_LDFLAGS=
ifdef I_MPI_ROOT
MPI_LDFLAGS+=-link_mpi=opt_mt
endif

# MPI Wrappers
WRAPPERS=I_MPI_CXX=$(MCXX) MPICH_CXX=$(MCXX) OMPI_CXX=$(MCXX)

# Extension name
ifeq ($(BSX),$(BSY))
EXT=$(BSX)bs.exe
else
EXT=$(BSX)x$(BSY)bs.exe
endif

# List of programs
PROGS=heat_seq.$(EXT)    \
    heat_ompss.$(EXT)    \
    heat_mpi.pure.$(EXT) \
    heat_mpi.omp.$(EXT)  \
    heat_mpi.task.$(EXT) \

ifdef INTEROPERABILITY_SRC
PROGS+=heat_mpi.interop.$(EXT)
endif

# Sources
SMP_SRC=src/common/misc.cpp src/smp/main.cpp
MPI_SRC=src/common/misc.cpp src/mpi/main.cpp

all: $(PROGS)

heat_seq.$(EXT): $(SMP_SRC) src/smp/solver_seq.cpp
	$(CXX) $(CPPFLAGS) $(CFLAGS) -o $@ $^ $(LDFLAGS)

heat_ompss.$(EXT): $(SMP_SRC) src/smp/solver_ompss.cpp
	$(MCXX) $(CPPFLAGS) $(MCCFLAGS) -o $@ $^ $(LDFLAGS)

heat_mpi.pure.$(EXT): $(MPI_SRC) src/mpi/solver_pure.cpp
	$(MPICXX) $(CPPFLAGS) $(CFLAGS) -o $@ $^ $(LDFLAGS) $(MPI_LDFLAGS)

heat_mpi.omp.$(EXT): $(MPI_SRC) src/mpi/solver_omp.cpp
	$(WRAPPERS) $(MPICXX) $(CPPFLAGS) $(MCCFLAGS) -o $@ $^ $(LDFLAGS) $(MPI_LDFLAGS)

heat_mpi.task.$(EXT): $(MPI_SRC) src/mpi/solver_task.cpp
	$(WRAPPERS) $(MPICXX) $(CPPFLAGS) $(MCCFLAGS) -o $@ $^ $(LDFLAGS) $(MPI_LDFLAGS)

heat_mpi.interop.$(EXT): $(MPI_SRC) src/mpi/solver_task.cpp interop/libmpiompss-interop.a
	$(WRAPPERS) $(MPICXX) -DINTEROPERABILITY $(CPPFLAGS) $(MCCFLAGS) -o $@ $^ $(LDFLAGS) $(MPI_LDFLAGS)

interop/libmpiompss-interop.a:
	$(MAKE) -C $(INTEROPERABILITY_SRC) -f Makefile.manual

check: all
	@./scripts/run-tests.sh $(PROGS)

clean:
	rm -f *.o *.exe

clean-all: clean
	$(MAKE) clean -C $(INTEROPERABILITY_SRC) -f Makefile.manual
