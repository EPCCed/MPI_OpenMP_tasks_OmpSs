##############################################################################
#
#  config.mk
#
#  Cray Archer
#
##############################################################################

ifeq ($(PE_ENV), CRAY)
  # Cray
  # msglevel_2 is cautions (rather tediously verbose)
  # msglevel_3 is warnings (ok)
  # Default C is c99

  CFLAGS_EXTRA=-hlist=m -h msglevel_3 -h stdc # -DNDEBUG -DKEEPHYDROONTARGET -DKEEPFIELDONTARGET #-h noomp
endif

ifeq ($(PE_ENV), INTEL)
  # Intel
  # -w2  gives errors and warnings
  # -w3  adds remarks (very verbose)

  CFLAGS_EXTRA= -w2 -strict-ansi -std=c99 -qopenmp
endif

ifeq ($(PE_ENV), GNU)
  # GNU
  CFLAGS_EXTRA= -g -Wall -pedantic -std=c99 -fopenmp -ftree-vectorizer-verbose=1 
endif

CFLAGS_EXTRA= -g -Wall -pedantic -std=c99 -fopenmp -ftree-vectorizer-verbose=1 

# CC=cc
# MPICC=cc
CC=mpicc
MPICC=mpicc
CFLAGS=-O2 $(CFLAGS_EXTRA)

LAUNCH_SERIAL_CMD=aprun -q -n 1
LAUNCH_MPI_CMD=aprun -q
LAUNCH_MPI_NP_SWITCH=-n
