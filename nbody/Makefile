# Compilers
CC=gcc
MPICC=mpicc
MCC=mcc

CPPFLAGS=
CPPFLAGS_SMP=
CFLAGS=
MCCFLAGS=
LDFLAGS=

# Nbody parameters
BIGO?=N2
BS?=2048

# Preprocessor flags
CPPFLAGS+=-Isrc -DBIGO=$(BIGO) -D_BIGO_$(BIGO) -DBLOCK_SIZE=$(BS)

# Compiler flags
CFLAGS+=-O3 -std=gnu11
MCCFLAGS+=--ompss-v2 $(CFLAGS) --Wn,-O3,-std=gnu11

# Linker flags
LDFLAGS+=-lrt -lm

# Compilation
PROGNAME=$(BIGO).$(BS).exe

# Directories
OBJDIR=obj/

all: makedir nbody_seq nbody_seq_blocking nbody_omp nbody_ompss nbody_mpi

nbody_seq: common/common.o plain/utils.o plain/solver_seq.o plain/main.o
	$(CC) $(CFLAGS) -o $@.$(PROGNAME) $(addprefix $(OBJDIR),$^) $(LDFLAGS)

nbody_seq_blocking: common/common.o blocking/common/utils.o blocking/smp/utils.o blocking/smp/solver_seq.o blocking/smp/main.o
	$(CC) $(CFLAGS) -o $@.$(PROGNAME) $(addprefix $(OBJDIR),$^) $(LDFLAGS)

nbody_omp: common/common.o plain/utils.o plain/solver_omp.o plain/main.o
	$(CC) -fopenmp $(CFLAGS) -o $@.$(PROGNAME) $(addprefix $(OBJDIR),$^) $(LDFLAGS)

nbody_ompss: common/common.o blocking/common/utils.o blocking/smp/utils.o blocking/smp/solver_ompss.o blocking/smp/main.o
	$(MCC) $(MCCFLAGS) -o $@.$(PROGNAME) $(addprefix $(OBJDIR),$^) $(LDFLAGS)

nbody_mpi: common/common.o blocking/common/utils.o blocking/mpi/utils.o blocking/mpi/solver_mpi.o blocking/mpi/main.o
	I_MPI_CC=$(MCC) MPICH_CC=$(MCC) OMPI_CC=$(MCC) $(MPICC) $(MCCFLAGS) -o $@.$(PROGNAME) $(addprefix $(OBJDIR),$^) $(LDFLAGS)

common/%.o: src/common/%.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c -o $(OBJDIR)$@ $^

plain/solver_omp.o: src/plain/solver_omp.c
	$(CC) -fopenmp $(CPPFLAGS) $(CFLAGS) -c -o $(OBJDIR)$@ $^

plain/%.o: src/plain/%.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c -o $(OBJDIR)$@ $^

blocking/smp/solver_ompss.o: src/blocking/smp/solver_ompss.c
	$(MCC) $(CPPFLAGS) $(MCCFLAGS) -c -o $(OBJDIR)$@ $^

blocking/smp/%.o: src/blocking/smp/%.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c -o $(OBJDIR)$@ $^

blocking/mpi/%.o: src/blocking/mpi/%.c
	I_MPI_CC=$(MCC) MPICH_CC=$(MCC) OMPI_CC=$(MCC) $(MPICC) $(CPPFLAGS) $(MCCFLAGS) -c -o $(OBJDIR)$@ $^

blocking/common/%.o: src/blocking/common/%.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c -o $(OBJDIR)$@ $^

makedir:
	mkdir -p obj/common obj/plain obj/blocking/common obj/blocking/smp obj/blocking/mpi

clean:
	rm -f *.o *.exe
	rm -rf obj/
