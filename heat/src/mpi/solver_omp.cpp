#include <mpi.h>
#include <algorithm>

#include "common/heat.hpp"

inline void solveBlock(block_t *matrix, int nbx, int nby, int bx, int by)
{
	block_t &targetBlock = matrix[bx*nby + by];
	const block_t &centerBlock = matrix[bx*nby + by];
	const block_t &topBlock    = matrix[(bx-1)*nby + by];
	const block_t &leftBlock   = matrix[bx*nby + (by-1)];
	const block_t &rightBlock  = matrix[bx*nby + (by+1)];
	const block_t &bottomBlock = matrix[(bx+1)*nby + by];
	
	double sum = 0.0;
	for (int x = 0; x < BSX; ++x) {
		const row_t &topRow = (x > 0) ? centerBlock[x-1] : topBlock[BSX-1];
		const row_t &bottomRow = (x < BSX-1) ? centerBlock[x+1] : bottomBlock[0];
		
		for (int y = 0; y < BSY; ++y) {
			double left = (y > 0) ? centerBlock[x][y-1] : leftBlock[x][BSY-1];
			double right = (y < BSY-1) ? centerBlock[x][y+1] : rightBlock[x][0];
			
			double value = 0.25 * (topRow[y] + bottomRow[y] + left + right);
			double diff = value - targetBlock[x][y];
			sum += diff * diff;
			targetBlock[x][y] = value;
		}
	}
}

inline void sendFirstComputeRow(block_t *matrix, int nbx, int nby, int rank, int rank_size)
{
	for (int by = 1; by < nby-1; ++by) {
		MPI_Send(&matrix[nby+by][0], BSY, MPI_DOUBLE, rank - 1, by, MPI_COMM_WORLD);
	}
}

inline void sendLastComputeRow(block_t *matrix, int nbx, int nby, int rank, int rank_size)
{
	for (int by = 1; by < nby-1; ++by) {
		MPI_Send(&matrix[(nbx-2)*nby + by][BSX-1], BSY, MPI_DOUBLE, rank + 1, by, MPI_COMM_WORLD);
	}
}

inline void receiveUpperBorder(block_t *matrix, int nbx, int nby, int rank, int rank_size)
{
	for (int by = 1; by < nby-1; ++by) {
		MPI_Recv(&matrix[by][BSX-1], BSY, MPI_DOUBLE, rank - 1, by, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
}

inline void receiveLowerBorder(block_t *matrix, int nbx, int nby, int rank, int rank_size)
{
	for (int by = 1; by < nby-1; ++by) {
		MPI_Recv(&matrix[(nbx-1)*nby + by][0], BSY, MPI_DOUBLE, rank + 1, by, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
}

inline void solveGaussSeidel(block_t *matrix, int nbx, int nby, int rank, int rank_size)
{
	if (rank != 0) {
		sendFirstComputeRow(matrix, nbx, nby, rank, rank_size);
		receiveUpperBorder(matrix, nbx, nby, rank, rank_size);
	}
	
	if (rank != rank_size - 1) {
		receiveLowerBorder(matrix, nbx, nby, rank, rank_size);
	}
	
	for (int bx = 1; bx < nbx-1; ++bx) {
		for (int by = 1; by < nby-1; ++by) {
			#pragma oss task label(gauss seidel)     \
					in(([nbx][nby]matrix)[bx-1][by]) \
					in(([nbx][nby]matrix)[bx][by-1]) \
					in(([nbx][nby]matrix)[bx][by+1]) \
					in(([nbx][nby]matrix)[bx+1][by]) \
					inout(([nbx][nby]matrix)[bx][by])
			solveBlock(matrix, nbx, nby, bx, by);
		}
	}
	#pragma oss taskwait
	
	if (rank != rank_size - 1) {
		sendLastComputeRow(matrix, nbx, nby, rank, rank_size);
	}
}

double solve(block_t *matrix, int rowBlocks, int colBlocks, int timesteps)
{
	int rank, rank_size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &rank_size);
	
	for (int t = 0; t < timesteps; ++t) {
		solveGaussSeidel(matrix, rowBlocks, colBlocks, rank, rank_size);
	}
	
	#pragma oss taskwait
	MPI_Barrier(MPI_COMM_WORLD);
	
	return 0.0;
}

