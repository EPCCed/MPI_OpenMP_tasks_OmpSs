#include "common.h"

#include <assert.h>
#include <getopt.h>
#include <ieee754.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <time.h>

void * nbody_alloc(size_t size)
{
	void *addr = mmap(NULL, size, PROT_READ|PROT_WRITE, MAP_PRIVATE|MAP_ANONYMOUS, -1, 0);
	assert(addr != MAP_FAILED);
	return addr;
}

nbody_conf_t nbody_get_conf(int argc, char **argv)
{
	nbody_conf_t conf;
	conf.domain_size_x = default_domain_size_x;
	conf.domain_size_y = default_domain_size_y;
	conf.domain_size_z = default_domain_size_z;
	conf.mass_maximum  = default_mass_maximum;
	conf.time_interval = default_time_interval;
	conf.seed          = default_seed;
	conf.name          = default_name;
	conf.num_particles = default_num_particles;
	conf.num_blocks    = conf.num_particles / BLOCK_SIZE;
	conf.timesteps     = default_timesteps;
	conf.save_result   = default_save_result;
	conf.check_result  = default_check_result;
	
	static struct option long_options[] = {
		{"particles",	required_argument,	0, 'p'},
		{"timesteps",	required_argument,	0, 't'},
		{"check",		no_argument,		0, 'c'},
		{"no-check",	no_argument,		0, 'C'},
		{"output",		no_argument,		0, 'o'},
		{"no-output",	no_argument,		0, 'O'},
		{"help",		no_argument,		0, 'h'},
		{0, 0, 0, 0}
	};
	
	int c;
	int index;
	while ((c = getopt_long(argc, argv, "hoOcCp:t:", long_options, &index)) != -1) {
		switch (c) {
			case 'h':
				fprintf(stderr, "Usage: %s [OPTION]...\n", argv[0]);
				fprintf(stderr, "Available options:\n");
				fprintf(stderr, "  -p, --particles=PARTICLES\t\tuse PARTICLES as the total number of particles (default: 16384)\n");
				fprintf(stderr, "  -t, --timesteps=TIMESTEPS\t\tuse TIMESTEPS as the number of timesteps (default: 10)\n");
				fprintf(stderr, "  -c, --check\t\t\t\tcheck the correctness of the result (disabled by default)\n");
				fprintf(stderr, "  -C, --no-check\t\t\tdo not check the correctness of the result\n");
				fprintf(stderr, "  -o, --output\t\t\t\tsave the computed particles to the default output file (disabled by default)\n");
				fprintf(stderr, "  -O, --no-output\t\t\tdo not save the computed particles to the default output file\n");
				fprintf(stderr, "  -h, --help\t\t\t\tdisplay this help and exit\n");
				exit(0);
			case 'o':
				conf.save_result = 1;
				break;
			case 'O':
				conf.save_result = 0;
				break;
			case 'c':
				conf.check_result = 1;
				break;
			case 'C':
				conf.check_result = 0;
				break;
			case 'p':
				conf.num_particles = atoi(optarg);
				break;
			case 't':
				conf.timesteps = atoi(optarg);
				break;
			case '?':
				exit(1);
			default:
				abort();
		}
	}
	return conf;
}

double nbody_compute_throughput(int num_particles, int timesteps, double elapsed_time)
{
	double interactions_per_timestep = 0;
#if defined(_BIGO_N2)
	interactions_per_timestep = (double)(num_particles)*(double)(num_particles);
#elif defined(_BIGO_NlogN)
	interactions_per_timestep = (double)(num_particles)*(double)(LOG2(num_particles));
#elif defined(_BIGO_N)
	interactions_per_timestep = (double)(num_particles);
#endif
	return (((interactions_per_timestep * (double)timesteps) / elapsed_time) / 1000000.0);
}

double get_time()
{
	struct timespec ts;
	clock_gettime(CLOCK_MONOTONIC, &ts);
	return (double)(ts.tv_sec) + (double)ts.tv_nsec * 1.0e-9;
}
