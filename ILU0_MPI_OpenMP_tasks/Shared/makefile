CC = icc
CC = gcc
CC = mpicc
CFLAGS = -Wall -g

AR = ar
ARFLAGS = ru

RL = ranlib

OBJS_CLOCK  = reloj.o 
OBJS_VECTOR = InputOutput.o IntVectors.o LongVectors.o FloatVectors.o DoubleVectors.o
OBJS_SPARSE = SparseVectors.o SparseMatrices.o SparseMatricesELL.o SparseSymmetric.o HarwellBoeing.o SparseHarwellBoeing.o 
OBJS_SPARSE_NEW = SparseVectorsNew.o SparseMatricesNew.o SparseSymmetricNew.o HarwellBoeing.o SparseHarwellBoeingNew.o 

OBJS = $(OBJS_CLOCK) $(OBJS_VECTOR) $(OBJS_SPARSE) $(OBJS_SPARSE_NEW)

default: libclock.a libvector.a libsparse.a libsparsenew.a

libshared.a : $(OBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RL) $(RLFLAGS) $@

libclock.a : $(OBJS_CLOCK)
	$(AR) $(ARFLAGS) $@ $?
	$(RL) $(RLFLAGS) $@

libvector.a : $(OBJS_VECTOR)
	$(AR) $(ARFLAGS) $@ $?
	$(RL) $(RLFLAGS) $@

libsparse.a : $(OBJS_SPARSE)
	$(AR) $(ARFLAGS) $@ $?
	$(RL) $(RLFLAGS) $@

libsparsenew.a : $(OBJS_SPARSE_NEW)
	$(AR) $(ARFLAGS) $@ $?
	$(RL) $(RLFLAGS) $@

.c.o:
	echo compilando
	$(CC) $(CFLAGS) -c $*.c

clean:
	rm *.o

clear:
	rm *.o *.a
