#include "plain/nbody.h"

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
#include <time.h>
#include <unistd.h>

void nbody_particle_init(const nbody_conf_t *conf, particles_t *part)
{
	const int num_particles = conf->num_particles;
	
	for (int i = 0; i < num_particles; i++) {
		part->position_x[i] = conf->domain_size_x * ((float)random() / ((float)RAND_MAX + 1.0));
		part->position_y[i] = conf->domain_size_y * ((float)random() / ((float)RAND_MAX + 1.0));
		part->position_z[i] = conf->domain_size_z * ((float)random() / ((float)RAND_MAX + 1.0));
		part->mass[i] = conf->mass_maximum * ((float)random() / ((float)RAND_MAX + 1.0));
		part->weight[i] = gravitational_constant * part->mass[i];
	}
}

int nbody_compare_particles(const particles_t *local, const particles_t *reference, int num_particles)
{
	double error = 0.0;
	int count = 0;
	for (int i = 0; i < num_particles; i++) {
		if ((local->position_x[i] != reference->position_x[i]) ||
		    (local->position_y[i] != reference->position_y[i]) ||
		    (local->position_z[i] != reference->position_z[i])) {
				error += fabs(((local->position_x[i] - reference->position_x[i])*100.0) / reference->position_x[i]) +
				         fabs(((local->position_y[i] - reference->position_y[i])*100.0) / reference->position_y[i]) +
				         fabs(((local->position_z[i] - reference->position_z[i])*100.0) / reference->position_z[i]);
				count++;
		}
	}
	
	double relative_error = error / (3.0 * count);
	if ((count * 100.0) / (num_particles) > 0.6 || relative_error >  0.000008) {
		printf("Relative error[%d]: %f\n", count, relative_error);
		return 0;
	}
	return 1;
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
	
	void *reference_addr = mmap(NULL, nbody->file.size, PROT_READ, MAP_SHARED, fd, 0);
	assert(reference_addr != MAP_FAILED);
	
	particles_t reference;
	nbody_link_particles(&reference, nbody->num_particles, reference_addr);
	
	if (nbody_compare_particles(&nbody->particles, &reference, nbody->num_particles)) {
		printf("Result validation: OK\n");
	}
	
	int err = munmap(reference_addr, nbody->file.size);
	assert(!err);
	
	err = close(fd);
	assert(!err);
}

void nbody_generate_particles(const nbody_conf_t *conf, nbody_file_t *file)
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
	
	void *addr = mmap(NULL, size, PROT_WRITE|PROT_READ, MAP_SHARED, fd, 0);
	assert(addr != MAP_FAILED);
	
	particles_t particles;
	nbody_link_particles(&particles, conf->num_particles, addr);
	
	// Initialize the particles
	nbody_particle_init(conf, &particles);
	
	err = munmap(addr, size);
	assert(!err);
	
	err = close(fd);
	assert(!err);
}

void nbody_load_particles(nbody_t *nbody, const nbody_conf_t *conf, nbody_file_t *file)
{
	char fname[1024];
	sprintf(fname, "%s.in", file->name);
	
	const int fd = open(fname, O_RDONLY, 0);
	assert(fd >= 0);
	
	void * addr = mmap(NULL, file->size, PROT_READ|PROT_WRITE, MAP_PRIVATE, fd, 0);
	assert(addr != MAP_FAILED);
	
	memcpy(nbody->particles._ptr, addr, nbody->particles._size);
	
	int err = munmap(addr, file->size);
	assert(!err);
	
	err = close(fd);
	assert(!err);
}

void nbody_setup_file(nbody_t *nbody, const nbody_conf_t *conf)
{
	nbody_file_t file;
	file.size = nbody->particles._size;
	sprintf(file.name, "%s-%s-%d-%d-%d", conf->name, TOSTRING(BIGO), conf->num_particles, conf->num_particles, conf->timesteps);
	nbody->file = file;
}

void nbody_link_particles(particles_t *particles, int num_particles, void *addr)
{
	const size_t size = num_particles * sizeof(float);
	particles->position_x = addr + (0 * size);
	particles->position_y = addr + (1 * size);
	particles->position_z = addr + (2 * size);
	particles->velocity_x = addr + (3 * size);
	particles->velocity_y = addr + (4 * size);
	particles->velocity_z = addr + (5 * size);
	particles->mass       = addr + (6 * size);
	particles->weight     = addr + (7 * size);
	particles->_ptr       = addr;
	particles->_size      = size * PARTICLE_MEMBERS;
}

void nbody_setup_particles(nbody_t *nbody, const nbody_conf_t *conf)
{
	void *addr = nbody_alloc(nbody->particles._size);
	nbody_link_particles(&nbody->particles, nbody->num_particles, addr);
	
	nbody_generate_particles(conf, &nbody->file);
	nbody_load_particles(nbody, conf, &nbody->file);
}

void nbody_setup_forces(nbody_t *nbody, const nbody_conf_t *conf)
{
	void *base_addr = nbody_alloc(nbody->forces._size);
	size_t array_size = conf->num_particles * sizeof(float);
	
	nbody->forces._ptr = base_addr;
	nbody->forces.x = base_addr + (0 * array_size);
	nbody->forces.y = base_addr + (1 * array_size);
	nbody->forces.z = base_addr + (2 * array_size);
}

nbody_t nbody_setup(const nbody_conf_t *conf)
{
	nbody_t nbody;
	nbody.num_particles = conf->num_particles;
	nbody.timesteps = conf->timesteps;
	
	// Compute the size of the structures
	nbody.particles._size = nbody.num_particles * PARTICLE_MEMBERS * sizeof(float);
	nbody.forces._size = nbody.num_particles * FORCE_MEMBERS * sizeof(float);
	
	// Set up the file
	nbody_setup_file(&nbody, conf);
	
	// Initialize the particles
	nbody_setup_particles(&nbody, conf);
	
	// Initialize the forces
	nbody_setup_forces(&nbody, conf);
	
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
	
	void *particles = mmap(NULL, size, PROT_WRITE|PROT_READ, MAP_SHARED, fd, 0);
	assert(particles != MAP_FAILED);
	
	memcpy(particles, nbody->particles._ptr, size);
	
	err = munmap(particles, size);
	assert(!err);
	
	err = close(fd);
	assert(!err);
}

void nbody_free(nbody_t *nbody)
{
	int err = munmap(nbody->particles._ptr, nbody->particles._size);
	err |= munmap(nbody->forces._ptr, nbody->forces._size);
	assert(!err);
}

