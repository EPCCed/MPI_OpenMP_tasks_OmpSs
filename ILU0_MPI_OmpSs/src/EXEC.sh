  
# For OPENMPI 1.8 or greater
#mpirun -np $1 --bind-to core --map-by node:pe=$2 --report-bindings ./Test  1.0e-2 5 10 $2 $3 $4

# For IMPI, MPICH
mpirun -np $1 ./Test  1.0e-2 5 10 $2 $3 $4


