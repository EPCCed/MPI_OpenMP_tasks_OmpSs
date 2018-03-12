#!/bin/bash

# MPI Execution configuration
nprocs=2
nthreadsxproc=2
totalthreads=$(($nprocs * $nthreadsxproc))

# Execution input
size=$((4*1024))
timesteps=100

if [ "$#" -lt 1 ]; then
	echo "Usage: $0 prog1 [prog2...]"
	exit 1
fi

programs=("$@")

export NANOS6=optimized

reffile=.heat_ref.ppm
tmpfile=.heat_tmp.ppm

echo ---------------------------------
echo TEST SUITE
echo ---------------------------------

total=0
failed=0

seq_prog=false
for prog in ${programs[*]}; do
	if [[ $prog == *"_seq"* ]] && [ -f ${prog} ]; then
		./${prog} -s $size -t $timesteps -o$reffile > /dev/null
		if [ $? -ne 0 ]; then
			echo "Sequential program has failed!"
			exit 1
		fi
		seq_prog=true
	fi
done

if [ "$seq_prog" == false ]; then
	echo "Sequential program not found!"
	exit 1
fi

for prog in ${programs[*]}; do
	# Check if it exists
	if [ ! -f ${prog} ]; then
		continue;
	fi
	total=$(($total + 1))

	# Execution
	if [[ $prog == *"_mpi.pure"* ]] || [[ $prog == *"_gaspi.pure"* ]]; then
		mpiexec.hydra -n $totalthreads ./${prog} -s $size -t $timesteps -o$tmpfile > /dev/null
	elif [[ $prog == *"_mpi"* ]]; then
		mpiexec.hydra -n $nprocs -bind-to hwthread:$nthreadsxproc ./${prog} -s $size -t $timesteps -o$tmpfile > /dev/null
	else
		./${prog} -s $size -t $timesteps -o$tmpfile > /dev/null
	fi

	# Check the return value of the program
	if [ $? -ne 0 ]; then
		failed=$(($failed + 1))
		echo $prog FAILED
		continue;
	fi

	# Result checking
	diff $reffile $tmpfile > /dev/null
	if [ $? -eq 0 ]; then
		echo $prog PASSED
	else
		failed=$(($failed + 1))
		echo $prog FAILED
	fi
done

echo ---------------------------------
if [ $failed -eq 0 ]; then
	echo SUMMARY: ALL TESTS PASSED
else
	echo SUMMARY: ${failed}/${total} TESTS FAILED
fi
echo ---------------------------------

rm -f $reffile
rm -f $tmpfile
