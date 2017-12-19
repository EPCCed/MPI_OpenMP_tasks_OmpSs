#!/bin/bash

# MPI Execution configuration
nprocs=2
nthreadsxproc=2

# Execution input
particles=$((16*1024))
timesteps=20

if [ "$#" -lt 1 ]; then
	echo "Usage: $0 prog1 [prog2...]"
	exit 1
fi

bigo=$1
bs=$2
programs=("$@")

export NANOS6=optimized
export NANOS6_SCHEDULER=fifo
export OMP_NUM_THREADS=$(nproc)

logfile=.nbody_log

echo ---------------------------------
echo TEST SUITE
echo ---------------------------------

total=0
failed=0

for prog in ${programs[*]}; do
	# Check if it exists
	if [ ! -f ${prog} ]; then
		    continue;
	fi

	# Execute with/without MPI
	if [[ $prog == *"_mpi"* ]]; then
		mpiexec.hydra -n $nprocs -bind-to hwthread:$nthreadsxproc ./${prog} -p $particles -t $timesteps -c &> $logfile
	else
		./${prog} -p $particles -t $timesteps -c &> $logfile
	fi

	# Result checking
	grep -Fxq 'Result validation: OK' $logfile
	if [ $? -eq 0 ]; then
		echo $prog PASSED
	else
		failed=$(($failed + 1))
		echo $prog FAILED
	fi
	total=$(($total + 1))
done

echo ---------------------------------
if [ $failed -eq 0 ]; then
	echo SUMMARY: ALL TESTS PASSED
else
	echo SUMMARY: ${failed}/${total} TESTS FAILED
fi
echo ---------------------------------

rm $logfile
