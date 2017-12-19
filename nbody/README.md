# N-body

## Introduction
An N-body simulation numerically approximates the evolution of a system of
bodies in which each body continuously interacts with every other body.  A
familiar example is an astrophysical simulation in which each body represents a
galaxy or an individual star, and the bodies attract each other through the
gravitational force.

N-body simulation arises in many other computational science problems as well.
For example, protein folding is studied using N-body simulation to calculate
electrostatic and van der Waals forces. Turbulent fluid flow simulation and
global illumination computation in computer graphics are other examples of
problems that use N-body simulation.

## Available versions and building instructions

The nbody application has several versions which are compiled in different 
binaries, by executing the `make` command. These versions can be blocking, 
when the particle space is divided into smaller blocks, or non-blocking, when 
it is not. They are:

  * **nbody_seq_plain**: Sequential version (non-blocking).
  * **nbody_omp_plain**: Parallel version using fork-join parallelism (non-blocking).
  * **nbody_seq**: Sequential version (blocking).
  * **nbody_ompss**: Parallel version using OmpSs tasks (blocking).
  * **nbody_mpi**: Parallel version using MPI (blocking).
  * **nbody_mpi_ompss**: Parallel version using OmpSs tasks + MPI (blocking). Communication tasks are serialized.
  * **nbody_mpi_ompss_interop**: Parallel version using OmpSs tasks + MPI + Interoperability library (blocking). *See building instructions, step 1*.

  The simplest way to compile this package is:

  1. `cd` to the directory containing the version of the benchmark
      that you want to build, or stay in N-Body root directory if
      you want to recursively build all the versions. 
      The N-Body + MPI + Interoperability library version is
      compiled only if the environment variable `INTEROPERABILITY_HOME`
      is set.

  2. Type `make` to compile the selected benchmark's version(s).
     Optionally, you can use a different block size when building
     the benchmark (by default 2048). Type `make BS=MY_BLOCK_SIZE`
     in order to change this value.

  3. In addition, you can type 'make check' to check the correctness
     of the built versions. By default, MPI versions run with 2 MPI
     processes and 2 hardware threads for each process, but you can
     change these parameters in 'scripts/run-tests.sh'.

## Execution instructions

The binaries accept several options. The most relevant options are the number 
of total particles with `-p` (default: 16384), and the number of timesteps with 
`-t` (default: 10). More options can be seen passing the `-h` option. An example 
of execution could be:

```
$ mpiexec -n 4 -bind-to hwthread:16 nbody_mpi_ompss.N2.1024.exe -t 100 -p 8192
```

in which the application will perform 100 timesteps in 4 MPI processes with 16 
hardware threads in each process (used by the OmpSs runtime). The total number 
of particles will be 8192, this means that each process will have 2048 
particles (2 blocks per process).

