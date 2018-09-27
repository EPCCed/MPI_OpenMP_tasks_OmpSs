include ../taudef  # include local settings

CFLAGS += -Wall
CFLAGS += -Wextra
CFLAGS += -Wshadow
CFLAGS += -O2
CFLAGS += -std=c99

ifeq ($(strip $(DEBUG_MODE)), 1)
  CFLAGS += -DDEBUG
endif

ifeq ($(strip $(DEBUG_MODE)), 2)
  CFLAGS += -O0 -ggdb
  CFLAGS += -DDEBUG
endif

ifeq ($(strip $(USE_MPI_MULTI_THREADED)), 1)
  CFLAGS += -DUSE_MPI_MULTI_THREADED
endif

################################################################################
INCLUDE_DIRS += $(NETCDF_INCDIR)
INCLUDE_DIRS += $(MPI_INCDIR)

LIBRARY_DIRS += $(NETCDF_LIBDIR)
LIBRARY_DIRS += $(MPI_LIBDIR)

LDFLAGS += $(addprefix -L,$(LIBRARY_DIRS))
CFLAGS += $(addprefix -I,$(INCLUDE_DIRS))
CFLAGS += -fopenmp

OBJ += exchange_matrix
OBJ += matrix
OBJ += linsolv
OBJ += linsys
OBJ += read_data
OBJ += setup_comm
OBJ += util

LIB += $(MPI_LIBNAME)
LIB += $(NETCDF_LIBNAME)
LIB += m

################################################################################

default: test_linsolv.exe

%linsolv.exe: %linsolv.o $(addsuffix .o, $(OBJ))
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) $(addprefix -l, $(LIB))

.INTERMEDIATE: linsolv.o

################################################################################

.PHONY: clean objclean

objclean:
	rm -f *.o

clean: objclean
	rm -f *.exe
