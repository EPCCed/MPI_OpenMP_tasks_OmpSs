#include "blocking/smp/nbody.h"

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


void nbody_generate_particles(const nbody_conf_t *conf, const nbody_file_t *file)
{
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
	
	const int size = file->size;
	assert(size % PAGE_SIZE == 0);
	
	int err = ftruncate(fd, size);
	assert(!err);
	
	particles_block_t * const particles = mmap(NULL, size, PROT_WRITE|PROT_READ, MAP_SHARED, fd, 0);
	
	for(int i = 0; i < conf->num_blocks; i++) {
		nbody_particle_init(conf, particles+i);
	}
	
	err = munmap(particles, size);
	assert(!err);
	
	err = close(fd);
	assert(!err);
}

void nbody_check(const nbody_t *nbody)
{
	char fname[1024];
	sprintf(fname, "%s.ref", nbody->file.name);
	if (access(fname, F_OK) != 0) {
		fprintf(stderr, "Warning: %s file does not exist. Skipping the check...\n", fname);
		return;
	}
	
	const int fd = open(fname, O_RDONLY, 0);
	assert(fd >= 0);
	
	particles_block_t *reference = mmap(NULL, nbody->file.size, PROT_READ, MAP_SHARED, fd, 0);
	assert(reference != MAP_FAILED);
	
	if (nbody_compare_particles(nbody->particles, reference, nbody->num_blocks)) {
		printf("Result validation: OK\n");
	} else {
		printf("Result validation: ERROR\n");
	}
	
	int err = munmap(reference, nbody->file.size);
	assert(!err);
	
	err = close(fd);
	assert(!err);
}

nbody_file_t nbody_setup_file(const nbody_conf_t *conf)
{
	nbody_file_t file;
	file.size = conf->num_blocks * sizeof(particles_block_t);
	
	sprintf(file.name, "%s-%s-%d-%d-%d", conf->name, TOSTRING(BIGO), conf->num_blocks * BLOCK_SIZE, BLOCK_SIZE, conf->timesteps);
	return file;
}

particles_block_t * nbody_load_particles(const nbody_conf_t *conf, const nbody_file_t *file)
{
	char fname[1024];
	sprintf(fname, "%s.in", file->name);
	
	const int fd = open(fname, O_RDONLY, 0);
	assert(fd >= 0);
	
	void * const ptr = mmap(NULL, file->size, PROT_READ|PROT_WRITE, MAP_PRIVATE, fd, 0);
	assert(ptr != MAP_FAILED);
	
	int err = close(fd);
	assert(!err);
	
	return ptr;
}

nbody_t nbody_setup(const nbody_conf_t *conf)
{
	nbody_file_t file = nbody_setup_file(conf);
	nbody_generate_particles(conf, &file);
	
	nbody_t nbody;
	nbody.particles = nbody_load_particles(conf, &file);
	assert(nbody.particles != NULL);
	
	nbody.forces = nbody_alloc(conf->num_blocks * sizeof(forces_block_t));
	assert(nbody.forces != NULL);
	
	nbody.num_blocks = conf->num_blocks;
	nbody.timesteps = conf->timesteps;
	nbody.file = file;
	
	return nbody;
}

void nbody_save_particles(const nbody_t *nbody)
{
	char fname[1024];
	sprintf(fname, "%s.out", nbody->file.name);
	
	const int fd = open(fname, O_RDWR|O_CREAT|O_TRUNC, S_IWUSR|S_IRUSR|S_IRGRP|S_IROTH);
	assert(fd >= 0);
	
	const int size = nbody->file.size;
	assert(size % PAGE_SIZE == 0);
	
	int err = ftruncate(fd, size);
	assert(!err);
	
	particles_block_t * const particles = mmap(NULL, size, PROT_WRITE|PROT_READ, MAP_SHARED, fd, 0);
	assert(particles != MAP_FAILED);
	
	memcpy(particles, nbody->particles, size);
	
	err = munmap(particles, size);
	assert(!err);
	
	err = close(fd);
	assert(!err);
}

void nbody_free(nbody_t *nbody)
{
	int err = munmap(nbody->particles, nbody->num_blocks * sizeof(particles_block_t));
	err |= munmap(nbody->forces, nbody->num_blocks * sizeof(forces_block_t));
	assert(!err);
}

