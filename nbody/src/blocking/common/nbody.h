#ifndef NBODY_H
#define NBODY_H

#include "common/common.h"

#include <stddef.h>
#include <stdlib.h>
#include <unistd.h>

// Block size definition
#ifndef BLOCK_SIZE
#define BLOCK_SIZE 2048
#endif

#define MIN_PARTICLES (4096 * BLOCK_SIZE / sizeof(particles_block_t))

// Solver structures
typedef struct {
	float position_x[BLOCK_SIZE]; /* m   */
	float position_y[BLOCK_SIZE]; /* m   */
	float position_z[BLOCK_SIZE]; /* m   */
	float velocity_x[BLOCK_SIZE]; /* m/s */
	float velocity_y[BLOCK_SIZE]; /* m/s */
	float velocity_z[BLOCK_SIZE]; /* m/s */
	float mass[BLOCK_SIZE];       /* kg  */
	float weight[BLOCK_SIZE];
} particles_block_t;

typedef struct {
	float x[BLOCK_SIZE]; /* x   */
	float y[BLOCK_SIZE]; /* y   */
	float z[BLOCK_SIZE]; /* z   */
} forces_block_t;

// Forward declaration
typedef struct nbody_file_t nbody_file_t;
typedef struct nbody_t nbody_t;

// Solver function
void nbody_solve(nbody_t *nbody, const int num_blocks, const int timesteps, const float time_interval);

// Auxiliary functions
nbody_t nbody_setup(const nbody_conf_t *conf);
void nbody_particle_init(const nbody_conf_t *conf, particles_block_t *part);
void nbody_stats(const nbody_t *nbody, const nbody_conf_t *conf, double time);
void nbody_save_particles(const nbody_t *nbody);
void nbody_free(nbody_t *nbody);
void nbody_check(const nbody_t *nbody);
int nbody_compare_particles(const particles_block_t *local, const particles_block_t *reference, int num_blocks);

#endif // NBODY_H

