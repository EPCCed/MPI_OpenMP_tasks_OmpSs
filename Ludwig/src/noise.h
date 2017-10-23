/*****************************************************************************
 *
 *  noise.h
 *
 *  Edinburgh Soft Matter and Statistical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *  Kevin Stratford (kevin@epcc.ed.ac.uk)
 *  (c) 2013 The University of Edinburgh
 *
 *****************************************************************************/

#ifndef NOISE_H
#define NOISE_H

#define NNOISE_MAX 10
#define NNOISE_STATE 4

#include "io_harness.h"

#include "targetDP.h"

typedef enum {NOISE_RHO = 0,
	      NOISE_PHI,
	      NOISE_QAB,
	      NOISE_END}
  noise_enum_t;

typedef struct noise_s noise_t;

__targetHost__ int noise_create(noise_t ** pobj);
__targetHost__ void noise_free(noise_t * obj);
__targetHost__ int noise_init(noise_t * obj, int master_seed);
__targetHost__ int noise_state_set(noise_t * obj, int index, unsigned int s[NNOISE_STATE]);
__targetHost__ int noise_state(noise_t * obj, int index, unsigned int s[NNOISE_STATE]);
__targetHost__  int noise_reap(noise_t * obj, int index, double * reap);
__targetHost__  int noise_reap_n(noise_t *obj, int index, int nmax, double * reap);
__targetHost__ int noise_uniform_double_reap(noise_t * obj, int index, double * reap);
__targetHost__ int noise_present(noise_t * obj, noise_enum_t type, int * present);
__targetHost__ int noise_present_set(noise_t * obj, noise_enum_t type, int present);
__targetHost__ int noise_init_io_info(noise_t * obj, int grid[3], int form_in, int form_out);

__targetHost__ unsigned int noise_uniform(unsigned int state[NNOISE_STATE]);

#endif
