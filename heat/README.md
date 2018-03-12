# Heat

## Introduction
The Heat simulation uses an iterative Gauss-Seidel method to solve the heat equation,
which is a parabolic partial differential equation that describes the distribution of
heat (or variation in temperature) in a given region over time.

The heat equation is of fundamental importance in a wide range of science fields. In
mathematics, it is the parabolic partial differential equation par excellence. In statistics,
it is related to the study of the Brownian motion. Also, the diffusion equation is a generic
version of the heat equation, and it is related to the study of chemical diffusion processes.

## Available versions and building instructions

The heat application has several versions which are compiled in different 
binaries, by executing the `make` command. They are:

  * **heat_seq**: Sequential version.
  * **heat_ompss**: Parallel version using OmpSs tasks.
  * **heat_mpi.pure**: Parallel version using MPI.
  * **heat_mpi.omp**: Parallel version using MPI + OmpSs. Fork-join parallelization.
  * **heat_mpi.task**: Parallel version using MPI + OmpSs tasks. Communication tasks are serialized.
  * **heat_mpi.interop**: Parallel version using MPI + OmpSs tasks + Interoperability library. *See building instructions, step 1*.


  The simplest way to compile this package is:

  1. Stay in Heat root directory to recursively build all the versions.
     The Heat MPI + OmpSs tasks + Interoperability library version is
     compiled only if the environment variable `INTEROPERABILITY_SRC`
     is set.

  2. Type `make` to compile the selected benchmark's version(s).
     Optionally, you can use a different block size in each dimension
     (BSX and BSY for vertical and horizontal dimensions, respectively)
     when building the benchmark (by default 1024). Type
     `make BSX=MY_BLOCK_SIZE_X BSY=MY_BLOCK_SIZE_Y` in order to change
     this value. If you want the same value in each dimension, type
     `make BSX=MY_BLOCK_SIZE`.

  3. In addition, you can type 'make check' to check the correctness
     of the built versions. By default, the pure MPI version runs with
     4 processes and the hybrid versions run with 2 MPI processes and 2
     hardware threads for each process. You can change these
     parameters in 'scripts/run-tests.sh'.

## Execution instructions

The binaries accept several options. The most relevant options are the size 
of the matrix with `-s` (default: 2048), and the number of timesteps with 
`-t` (default: 100). More options can be seen passing the `-h` option. An example 
of execution could be:

```
$ mpiexec -n 4 -bind-to hwthread:16 heat_mpi.task.1024bs.exe -t 150 -s 8192
```

in which the application will perform 150 timesteps in 4 MPI processes with 16 
hardware threads in each process (used by the OmpSs runtime). The size of the
matrix in each dimension will be 8192 (8192^2 elements in total), this means
that each process will have 2048 * 8192 elements (16 blocks per process).

