#include <mpi.h>
#include <cassert>
#include <iostream>

#include "common/heat.hpp"

#ifdef _OMPSS_2
#include <nanos6/debug.h>
#endif

#ifdef INTEROPERABILITY
#define DESIRED_THREAD_LEVEL (MPI_THREAD_MULTIPLE+1)
#else
#define DESIRED_THREAD_LEVEL (MPI_THREAD_MULTIPLE)
#endif

void generateImage(const HeatConfiguration &conf, int rowBlocks, int colBlocks, int rowBlocksPerRank);

int main(int argc, char **argv)
{
	int provided;
	MPI_Init_thread(&argc, &argv, DESIRED_THREAD_LEVEL, &provided);
	assert(provided == DESIRED_THREAD_LEVEL);
	
	int rank, rank_size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &rank_size);
	
	HeatConfiguration conf = readConfiguration(argc, argv);
	refineConfiguration(conf, rank_size * BSX, BSY);
	if (!rank) printConfiguration(conf);
	
	conf.rowBlocks = conf.rows / BSX;
	conf.colBlocks = conf.cols / BSY;
	int rowBlocks = conf.rowBlocks + 2;
	int colBlocks = conf.colBlocks + 2;
	int rowBlocksPerRank = conf.rowBlocks / rank_size + 2;
	
	int err = initialize(conf, rowBlocksPerRank, colBlocks, (rowBlocksPerRank - 2) * rank);
	assert(!err);
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	// Solve the problem
	double start = get_time();
	solve(conf.matrix, rowBlocksPerRank, colBlocks, conf.timesteps);
	double end = get_time();
	
	if (!rank) {
		long totalElements = (long)conf.rows * (long)conf.cols;
		double performance = totalElements * (long)conf.timesteps;
		performance = performance / (end - start);
		performance = performance / 1000000.0;
		
#ifdef _OMPSS_2
		int threads = nanos_get_num_cpus();
#else
		int threads = 1;
#endif
		
		fprintf(stdout, "rows, %d, cols, %d, rows_per_rank, %d, total, %ld, total_per_rank, %ld, bs, %d"
				" ,ranks, %d, threads, %d, timesteps, %d, time, %f, performance, %f\n",
				conf.rows, conf.cols, conf.rows / rank_size, totalElements, totalElements / rank_size,
				BSX, rank_size, threads, conf.timesteps, end - start, performance);
	}
	
	if (conf.generateImage) {
		generateImage(conf, rowBlocks, colBlocks, rowBlocksPerRank);
	}
	
	err = finalize(conf);
	assert(!err);
	
	MPI_Finalize();
	
	return 0;
}

void generateImage(const HeatConfiguration &conf, int rowBlocks, int colBlocks, int rowBlocksPerRank)
{
	int rank, rank_size;
	block_t *auxMatrix = nullptr;
	
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &rank_size);
	
	if (!rank) {
		auxMatrix = (block_t *) malloc(rowBlocks * colBlocks * sizeof(block_t));
		if (auxMatrix == NULL) {
			fprintf(stderr, "Memory cannot be allocated!\n");
			exit(1);
		}
		
		initializeMatrix(conf, auxMatrix, rowBlocks, colBlocks);
	}
	
	int count = (rowBlocksPerRank - 2) * colBlocks * BSX * BSY;
	MPI_Gather(
		&conf.matrix[colBlocks], count, MPI_DOUBLE,
		&auxMatrix[colBlocks], count, MPI_DOUBLE,
		0, MPI_COMM_WORLD
	);
	
	if (!rank) {
		int err = writeImage(conf.imageFileName, auxMatrix, rowBlocks, colBlocks);
		assert(!err);
		
		free(auxMatrix);
	}
}
