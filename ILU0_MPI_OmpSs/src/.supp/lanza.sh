#!/bin/bash

#SBATCH --job-name=ILUPACK
#SBATCH -D .
#SBATCH --exclusive
#SBATCH --output=ILUPACK_%j.out
#SBATCH --error=ILUPACK_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=1
#SBATCH --time=00:30:00

cd /gpfs/scratch/bsc19/bsc19844/intertwine/ILU0_mpi_ompss_Nanos6/src

module purge
module load intel
module load ompss
module load openmpi
#module load EXTRAE/latest
procesos_ompss=8
export NX_ARGS="--summary --verbose --smp-workers=$procesos_ompss" 

mpirun --bind-to core --map-by node:pe=8 ./Test  1.0e-2 5 10 1 8 8 ../../../A100.rsa 

