#include <algorithm>

#include "common/heat.hpp"


inline double solveBlock(block_t *matrix, int nbx, int nby, int bx, int by)
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
		
		//#pragma simd reduction(+:sum) vectorlength(64)
		for (int y = 0; y < BSY; ++y) {
			double left = (y > 0) ? centerBlock[x][y-1] : leftBlock[x][BSY-1];
			double right = (y < BSY-1) ? centerBlock[x][y+1] : rightBlock[x][0];
			
			double value = 0.25 * (topRow[y] + bottomRow[y] + left + right);
			double diff = value - targetBlock[x][y];
			sum += diff * diff;
			targetBlock[x][y] = value;
		}
	}
	
	return sum;
}

inline double gaussSeidelSolver(block_t *matrix, int nbx, int nby)
{
	double unew, diff, sum = 0.0;
	
	for (int bx = 1; bx < nbx-1; ++bx) {
		for (int by = 1; by < nby-1; ++by) {
			sum += solveBlock(matrix, nbx, nby, bx, by);
		}
	}
	
	return sum;
}

double solve(block_t *matrix, int rowBlocks, int colBlocks, int timesteps)
{
	double residual = 0.0;
	
	for (int t = 0; t < timesteps; ++t) {
		residual = gaussSeidelSolver(matrix, rowBlocks, colBlocks);
	}
	
	return residual;
}

