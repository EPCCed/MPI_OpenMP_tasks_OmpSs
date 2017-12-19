# Compilers
CC=gcc
MCC=mcc
MPICC=mpicc
MPI_WRAPPERS=I_MPI_CC=$(MCC) MPICH_CC=$(MCC) OMPI_CC=$(MCC)

CPPFLAGS=
CFLAGS=
MCCFLAGS=
LDFLAGS=

# Nbody parameters
BIGO?=N2
BS?=2048

# Interoperability flags
INTEROP_FLAGS=$(INTEROPERABILITY_HOME)/lib/libmpiompss-interop.a

# Preprocessor flags
CPPFLAGS+=-Isrc -DBIGO=$(BIGO) -D_BIGO_$(BIGO) -DBLOCK_SIZE=$(BS)

# Compiler flags
CFLAGS+=-O3 -std=gnu11
MCCFLAGS+=--ompss-2 $(CFLAGS) --Wn,-O3,-std=gnu11

# Linker flags
LDFLAGS+=-lrt -lm -lstdc++

SEQUENTIAL_PLAIN_SOURCES=  \
	src/common/common.c    \
	src/plain/utils.c      \
	src/plain/solver_seq.c \
	src/plain/main.c

SEQUENTIAL_SOURCES=                    \
	src/common/common.c                \
	src/blocking/common/common_utils.c \
	src/blocking/smp/utils.c           \
	src/blocking/smp/solver_seq.c      \
	src/blocking/smp/main.c

OPENMP_PLAIN_SOURCES=      \
	src/common/common.c    \
	src/plain/utils.c      \
	src/plain/solver_omp.c \
	src/plain/main.c

OMPSS_SOURCES=                         \
	src/common/common.c                \
	src/blocking/common/common_utils.c \
	src/blocking/smp/utils.c           \
	src/blocking/smp/solver_ompss.c    \
	src/blocking/smp/main.c

MPI_SOURCES=                           \
	src/common/common.c                \
	src/blocking/common/common_utils.c \
	src/blocking/mpi/utils.c           \
	src/blocking/mpi/solver_mpi.c      \
	src/blocking/mpi/main.c

MPI_OMPSS_SOURCES=                      \
	src/common/common.c                 \
	src/blocking/common/common_utils.c  \
	src/blocking/mpi/utils.c            \
	src/blocking/mpi/solver_mpi_ompss.c \
	src/blocking/mpi/main.c

PROGS= \
	nbody_seq_plain.$(BIGO).exe       \
	nbody_omp_plain.$(BIGO).exe       \
	nbody_seq.$(BIGO).$(BS).exe       \
	nbody_ompss.$(BIGO).$(BS).exe     \
	nbody_mpi.$(BIGO).$(BS).exe       \
	nbody_mpi_ompss.$(BIGO).$(BS).exe

ifdef INTEROPERABILITY_HOME
PROGS+=nbody_mpi_ompss_interop.$(BIGO).$(BS).exe
endif

all: $(PROGS)

nbody_seq_plain.$(BIGO).exe: $(SEQUENTIAL_PLAIN_SOURCES)
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^ $(LDFLAGS)

nbody_omp_plain.$(BIGO).exe: $(OPENMP_PLAIN_SOURCES)
	$(CC) $(CPPFLAGS) $(CFLAGS) -fopenmp -o $@ $^ $(LDFLAGS)

nbody_seq.$(BIGO).$(BS).exe: $(SEQUENTIAL_SOURCES)
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^ $(LDFLAGS)

nbody_ompss.$(BIGO).$(BS).exe: $(OMPSS_SOURCES)
	$(MCC) $(CPPFLAGS) $(MCCFLAGS) -o $@ $^ $(LDFLAGS)

nbody_mpi.$(BIGO).$(BS).exe: $(MPI_SOURCES)
	$(MPICC) $(CPPFLAGS) $(CFLAGS) -o $@ $^ $(LDFLAGS)

nbody_mpi_ompss.$(BIGO).$(BS).exe: $(MPI_OMPSS_SOURCES)
	$(MPI_WRAPPERS) $(MPICC) $(CPPFLAGS) $(MCCFLAGS) -o $@ $^ $(LDFLAGS)

nbody_mpi_ompss_interop.$(BIGO).$(BS).exe: $(MPI_OMPSS_SOURCES)
	$(MPI_WRAPPERS) $(MPICC) -DINTEROPERABILITY $(CPPFLAGS) $(MCCFLAGS) -o $@ $^ $(LDFLAGS) $(INTEROP_FLAGS)

check: all
	@./scripts/run-tests.sh $(PROGS)

clean:
	rm -f *.o *.exe
