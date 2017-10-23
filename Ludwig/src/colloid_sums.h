/*****************************************************************************
 *
 *  colloid_sums.h
 *
 *  $Id: colloid_sums.h 2479 2014-09-25 11:24:42Z stratford $
 *
 *  Edinburgh Soft Matter and Statistical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *  Kevin Stratford (kevin@epcc.ed.ac.uk)
 *  (c) 2010 The University of Edinburgh
 *
 *****************************************************************************/

#ifndef COLLOID_SUMS_H
#define COLLOID_SUMS_H

#include "colloids.h"

typedef struct colloid_sum_s colloid_sum_t;

typedef enum colloid_sum_enum_type {
  COLLOID_SUM_STRUCTURE = 0,
  COLLOID_SUM_DYNAMICS = 1,
  COLLOID_SUM_ACTIVE = 2,
  COLLOID_SUM_SUBGRID = 2,
  COLLOID_SUM_CONSERVATION = 3,
  COLLOID_SUM_MAX = 4} colloid_sum_enum_t;

/* Note that ACTIVE and SUBGRID are indeed the same */

int colloid_sums_create(colloids_info_t * cinfo, colloid_sum_t ** psum);
void colloid_sums_free(colloid_sum_t * sum);
int colloid_sums_halo(colloids_info_t * cinfo, colloid_sum_enum_t type);
int colloid_sums_1d(colloid_sum_t * sum, int dim, colloid_sum_enum_t type);

#endif
