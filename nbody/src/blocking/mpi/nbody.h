#ifndef NBODY_MPI_H
#define NBODY_MPI_H

#include "blocking/common/nbody.h"

// Application structures
struct nbody_file_t {
	size_t total_size;
	size_t size;
	size_t offset;
	char name[1000];
};

struct nbody_t {
	particles_block_t *local;
	particles_block_t *remote1;
	particles_block_t *remote2;
	forces_block_t *forces;
	int num_blocks;
	int timesteps;
	nbody_file_t file;
};

#endif // NBODY_MPI_H

