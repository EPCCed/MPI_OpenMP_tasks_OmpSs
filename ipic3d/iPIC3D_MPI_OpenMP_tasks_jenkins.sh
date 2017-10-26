#This release was prepared by Dana Akhmetova <danaak@kth.se>/<danieka@gmail.com> on behalf of the INTERTWinE European Exascale Project <http://www.intertwine-project.eu>
#!/bin/bash

d=$(pwd)

cd ipic3d/MPI_OpenMP_tasks
ipic3d_dir=$(pwd)
cd $d

rm -rf build-ipic3d
mkdir build-ipic3d
mkdir build-ipic3d/data
cd build-ipic3d

cmake $ipic3d_dir
make
ctest -V

#mpirun -n 6 ./iPic3D $ipic3d_dir/inputfiles/testMagnetosphere2Dsmall.inp
#diff data/ConservedQuantities.txt $ipic3d_dir/testfiles/ConservedQuantities.txt-testMagnetosphere2Dsmall
#bash -x $ipic3d_dir/scripts/compareconserved
