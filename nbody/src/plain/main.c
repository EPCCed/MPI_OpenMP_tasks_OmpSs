#include "plain/nbody.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>


int main(int argc, char** argv)
{
	nbody_conf_t conf = nbody_get_conf(argc, argv);
	
	conf.num_blocks = 1;
	conf.num_particles = ROUNDUP(conf.num_particles, MIN_PARTICLES);
	assert(conf.num_particles >= MIN_PARTICLES);
	assert(conf.timesteps > 0);
	
	nbody_t nbody = nbody_setup(&conf);
	
	double start = get_time();
	nbody_solve(&nbody, conf.timesteps, conf.time_interval);
	double end = get_time();
	
	nbody_stats(&nbody, &conf, end - start);
	
	if (conf.save_result) nbody_save_particles(&nbody);
	if (conf.check_result) nbody_check(&nbody);
	nbody_free(&nbody);
	return 0;
}
