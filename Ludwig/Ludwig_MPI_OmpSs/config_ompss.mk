##############################################################################
#
#  config.mk
#
#  Cray Archer
#
##############################################################################


ifeq ($(PE_ENV), GNU)
  # GNU
  CFLAGS_EXTRA= -Wall -std=c99 
endif


#CC=mpicc	
#MPICC=mpicc

CFLAGS_EXTRA= -Wall -std=c99 	
CC=/lustre/home/z04/lcebaman/ompss/bin/imcc 
MPICC=${CC}
CFLAGS=-O2 $(CFLAGS_EXTRA)

