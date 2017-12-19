#include "blocking/mpi/nbody.h"

#include <assert.h>
#include <ctype.h>
#include <fcntl.h>
#include <getopt.h>
#include <ieee754.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>

#include <mpi.h>


void nbody_generate_particles(const nbody_conf_t *conf, const nbody_file_t *file)
{
	int rank_size;
	MPI_Comm_size(MPI_COMM_WORLD, &rank_size);
	
	char fname[1024];
	sprintf(fname, "%s.in", file->name);
	
	if (!access(fname, F_OK)) {
		return;
	}
	
	struct stat st = {0};
	if (stat("data", &st) == -1) {
		mkdir("data", 0755);
	}
	
	const int fd = open(fname, O_RDWR|O_CREAT|O_TRUNC, S_IRUSR|S_IRGRP|S_IROTH);
	assert(fd >= 0);
	
	const int total_size = file->total_size;
	assert(total_size % PAGE_SIZE == 0);
	
	int err = ftruncate(fd, total_size);
	assert(!err);
	
	particles_block_t * const particles = mmap(NULL, total_size, PROT_WRITE|PROT_READ, MAP_SHARED, fd, 0);
	
	for(int i = 0; i < conf->num_blocks * rank_size; i++) {
		nbody_particle_init(conf, particles+i);
	}
	
	err = munmap(particles, total_size);
	assert(!err);
	
	err = close(fd);
	assert(!err);
}

void nbody_check(const nbody_t *nbody)
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	char fname[1024];
	sprintf(fname, "%s.ref", nbody->file.name);
	if (access(fname, F_OK) != 0) {
		if (!rank) fprintf(stderr, "Warning: %s file does not exist. Skipping the check...\n", fname);
		return;
	}
	
	const int fd = open(fname, O_RDONLY, 0);
	assert(fd >= 0);
	
	particles_block_t *reference = mmap(NULL, nbody->file.size, PROT_READ, MAP_SHARED, fd, nbody->file.offset);
	assert(reference != MAP_FAILED);
	
	int correct, correct_chunk;
	correct_chunk = nbody_compare_particles(nbody->local, reference, nbody->num_blocks);
	
	MPI_Reduce(&correct_chunk, &correct, 1, MPI_INT, MPI_LAND, 0, MPI_COMM_WORLD);
	
	if (!rank) {
		if (correct) {
			printf("Result validation: OK\n");
		} else {
			printf("Result validation: ERROR\n");
		}
	}
	
	int err = munmap(reference, nbody->file.size);
	assert(!err);
	
	err = close(fd);
	assert(!err);
}

nbody_file_t nbody_setup_file(const nbody_conf_t *conf)
{
	int rank, rank_size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &rank_size);
	
	nbody_file_t file;
	file.size = conf->num_blocks * sizeof(particles_block_t);
	file.total_size = rank_size * file.size;
	file.offset = rank * file.size;
	
	sprintf(file.name, "%s-%s-%d-%d-%d", conf->name, TOSTRING(BIGO), rank_size * conf->num_blocks * BLOCK_SIZE, BLOCK_SIZE, conf->timesteps);
	return file;
}

particles_block_t * nbody_load_particles(const nbody_conf_t *conf, const nbody_file_t *file)
{
	char fname[1024];
	sprintf(fname, "%s.in", file->name);
	
	const int fd = open(fname, O_RDONLY, 0);
	assert(fd >= 0);
	
	void * const ptr = mmap(NULL, file->size, PROT_READ|PROT_WRITE, MAP_PRIVATE, fd, file->offset);
	assert(ptr != MAP_FAILED);
	
	int err = close(fd);
	assert(!err);
	
	return ptr;
}

nbody_t nbody_setup(const nbody_conf_t *conf)
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	nbody_file_t file = nbody_setup_file(conf);
	
	if (!rank) nbody_generate_particles(conf, &file);
	MPI_Barrier(MPI_COMM_WORLD);
	
	nbody_t nbody;
	nbody.local = nbody_load_particles(conf, &file);
	assert(nbody.local != NULL);
	
	nbody.remote1 = nbody_alloc(conf->num_blocks * sizeof(particles_block_t));
	nbody.remote2 = nbody_alloc(conf->num_blocks * sizeof(particles_block_t));
	assert(nbody.remote1 != NULL);
	assert(nbody.remote2 != NULL);
	
	nbody.forces = nbody_alloc(conf->num_blocks * sizeof(forces_block_t));
	assert(nbody.forces != NULL);
	
	nbody.num_blocks = conf->num_blocks;
	nbody.timesteps = conf->timesteps;
	nbody.file = file;
	
	return nbody;
}

void nbody_save_particles(const nbody_t *nbody)
{
	int rank, rank_size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &rank_size);
	
	char fname[1024];
	sprintf(fname, "%s.out", nbody->file.name);
	const int mode = MPI_MODE_CREATE | MPI_MODE_WRONLY;
	
	MPI_File outfile;
	int err = MPI_File_open(MPI_COMM_WORLD, fname, mode, MPI_INFO_NULL, &outfile);
	assert(err == MPI_SUCCESS);
	
	MPI_File_set_view(outfile, nbody->file.offset, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);
	assert(err == MPI_SUCCESS);
	
	MPI_File_write(outfile, nbody->local, nbody->file.size, MPI_BYTE, MPI_STATUS_IGNORE);
	assert(err == MPI_SUCCESS);
	
	MPI_File_close(&outfile);
	assert(err == MPI_SUCCESS);
	
	MPI_Barrier(MPI_COMM_WORLD);
}

void nbody_free(nbody_t *nbody)
{
	const int particles_size = nbody->num_blocks * sizeof(particles_block_t);
	const int forces_size = nbody->num_blocks * sizeof(forces_block_t);
	
	int err = munmap(nbody->local, particles_size);
	err |= munmap(nbody->remote1, particles_size);
	err |= munmap(nbody->remote2, particles_size);
	err |= munmap(nbody->forces, forces_size);
	assert(!err);
}

