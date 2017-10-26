#!/bin/bash

MPIRUN=${mympirun:-mpirun}
TL=${taskloop:-0}

TEST_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR=$(dirname "${TEST_DIR}")
SRC_DIR=${ROOT_DIR}/src
DATA_DIR=${ROOT_DIR}/data

prog=test_linsolv.exe
tcase=NACA64A010
proclist="1 4"

declare -i status=0

function cleanup()
{
  rm -f ${TEST_DIR}/*.log.test
  cd ${SRC_DIR} && make -s clean
}

################################################################################
cleanup

# build linsolv
cd ${SRC_DIR}
echo "### Build linsolv ###"
CFLAGS+="-DUSE_TASKLOOP=${TL}" make -s DEBUG_MODE=0 ${prog}
  
if [ $? -eq 0 ]; then
  echo " Build successful"
else
  echo " Build failed!"
  status=$(($status + 1))
fi

# run linsolv
cd ${DATA_DIR}/${tcase}
export OMP_NUM_THREADS=4
for np in $proclist
do 
  log=${TEST_DIR}/test_linsolv_${tcase}_${np}.log.test

  echo "### Run ${tcase} with ${np} procs ###"
  ${MPIRUN} -np ${np} ${SRC_DIR}/${prog} > ${log}

  if [ $? -eq 0 ]; then
    echo " Execution successful"
  else
    echo " Execution failed!"
    status=$(($status + 1))
  fi
done

# compare results
cd ${TEST_DIR}
for np in $proclist
do
  reflog=test_linsolv_${tcase}_${np}.log.ref
  testlog=${reflog/.ref/.test}
  echo "### Comparing results for ${np} procs ###"
  if [ ! -e ${reflog} ]; then
    echo " ${reflog} is missing!"
    status=$(($status + 1))
    continue
  fi
  if [ ! -e ${testlog} ]; then
    echo " ${testlog} is missing!"
    status=$(($status + 1))
    continue
  fi

  diffoutput=$(diff -I "\[0\] Active .* version: " ${reflog} ${testlog})
  if [ $? -eq 0 ]; then
    echo " No differences"
  else
    echo "diff ${reflog} ${testlog}"
    echo "${diffoutput}"
    status=$(($status + 1))
  fi
done

if [ $status -eq 0 ]; then
  cleanup
fi

exit $status