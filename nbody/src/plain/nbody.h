#ifndef NBODY_H
#define NBODY_H

#include "common/common.h"

#include <stddef.h>
#include <stdlib.h>
#include <unistd.h>

#define MIN_PARTICLES 4096
#define PARTICLE_MEMBERS 8
#define FORCE_MEMBERS 3

// Solver structures
typedef struct {
	float *position_x; /* m   */
	float *position_y; /* m   */
	float *position_z; /* m   */
	float *velocity_x; /* m/s */
	float *velocity_y; /* m/s */
	float *velocity_z; /* m/s */
	float *mass;       /* kg  */
	float *weight;
	void *_ptr;
	size_t _size;
} particles_t;

typedef struct {
	float *x; /* x   */
	float *y; /* y   */
	float *z; /* z   */
	void *_ptr;
	size_t _size;
} forces_t;

// Application structures
typedef struct {
	size_t size;
	char name[1000];
} nbody_file_t;

typedef struct {
	particles_t particles;
	forces_t forces;
	int num_particles;
	int timesteps;
	nbody_file_t file;
} nbody_t;

// Solver function
void nbody_solve(nbody_t *nbody, const int timesteps, const float time_interval);

// Auxiliary functions
nbody_t nbody_setup(const nbody_conf_t *conf);
void nbody_setup_particles(nbody_t *nbody, const nbody_conf_t *conf);
void nbody_setup_forces(nbody_t *nbody, const nbody_conf_t *conf);
void nbody_link_particles(particles_t *reference, int num_particles, void *addr);
void nbody_particle_init(const nbody_conf_t *conf, particles_t *part);
void nbody_stats(const nbody_t *nbody, const nbody_conf_t *conf, double time);
void nbody_save_particles(const nbody_t *nbody);
void nbody_free(nbody_t *nbody);
void nbody_check(const nbody_t *nbody);
int nbody_compare_particles(const particles_t *local, const particles_t* reference, int num_particles);

#endif // NBODY_H

