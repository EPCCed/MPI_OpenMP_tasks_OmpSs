#include "plain/nbody.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


static void calculate_forces_N2(forces_t *forces, const particles_t *particles, const int num_particles);
static void update_particles(particles_t *particles, forces_t *forces, const int num_particles, const float time_interval);

void nbody_solve(nbody_t *nbody, const int timesteps, const float time_interval)
{
	assert(nbody != NULL);
	
	int num_particles = nbody->num_particles;
	particles_t *particles = &nbody->particles;
	forces_t *forces = &nbody->forces;
	
	for (int t = 0; t < timesteps; t++) {
		calculate_forces(forces, particles, num_particles);
		update_particles(particles, forces, num_particles, time_interval);
	}
}

void calculate_forces_N2(forces_t *forces, const particles_t *particles, const int num_particles)
{
	float *x = forces->x;
	float *y = forces->y;
	float *z = forces->z;
	
	const float *pos_x = particles->position_x;
	const float *pos_y = particles->position_y;
	const float *pos_z = particles->position_z;
	const float *mass  = particles->mass;
	
	for (int i = 0; i < num_particles; i++) {
		float fx = x[i], fy = y[i], fz = z[i];
		for (int j = 0; j < num_particles; j++) {
			const float diff_x = pos_x[j] - pos_x[i];
			const float diff_y = pos_y[j] - pos_y[i];
			const float diff_z = pos_z[j] - pos_z[i];
			
			const float distance_squared = diff_x * diff_x + diff_y * diff_y + diff_z * diff_z;
			const float distance = sqrtf(distance_squared);
			
			float force = 0.0f;
			if (distance_squared != 0.0f) {
				force = (mass[i] / (distance_squared * distance)) * (mass[j] * gravitational_constant);
			}
			fx += force * diff_x;
			fy += force * diff_y;
			fz += force * diff_z;
		}
		x[i] = fx;
		y[i] = fy;
		z[i] = fz;
	}
}

void update_particles(particles_t *particles, forces_t *forces, const int num_particles, const float time_interval)
{
	for (int e = 0; e < num_particles; e++) {
		const float mass       = particles->mass[e];
		const float velocity_x = particles->velocity_x[e];
		const float velocity_y = particles->velocity_y[e];
		const float velocity_z = particles->velocity_z[e];
		const float position_x = particles->position_x[e];
		const float position_y = particles->position_y[e];
		const float position_z = particles->position_z[e];
		
		const float time_by_mass       = time_interval / mass;
		const float half_time_interval = 0.5f * time_interval;
		
		const float velocity_change_x = forces->x[e] * time_by_mass;
		const float velocity_change_y = forces->y[e] * time_by_mass;
		const float velocity_change_z = forces->z[e] * time_by_mass;
		const float position_change_x = velocity_x + velocity_change_x * half_time_interval;
		const float position_change_y = velocity_y + velocity_change_y * half_time_interval;
		const float position_change_z = velocity_z + velocity_change_z * half_time_interval;
		
		particles->velocity_x[e] = velocity_x + velocity_change_x;
		particles->velocity_y[e] = velocity_y + velocity_change_y;
		particles->velocity_z[e] = velocity_z + velocity_change_z;
		particles->position_x[e] = position_x + position_change_x;
		particles->position_y[e] = position_y + position_change_y;
		particles->position_z[e] = position_z + position_change_z;
	}
	
	memset(forces->_ptr, 0, forces->_size);
}

void nbody_stats(const nbody_t *nbody, const nbody_conf_t *conf, double time)
{
	printf("bigo, %s, timesteps, %d, total_particles, %d, time, %g, performance, %g\n",
			TOSTRING(BIGO), nbody->timesteps, nbody->num_particles, time,
			nbody_compute_throughput(nbody->num_particles, nbody->timesteps, time)
	);
}
