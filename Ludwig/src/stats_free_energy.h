/****************************************************************************
 *
 *  stats_free_energy.h
 *
 *  $Id: stats_free_energy.h 2479 2014-09-25 11:24:42Z stratford $
 *
 *  Edinburgh Soft Matter and Statistical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *  Kevin Stratford (kevin@epcc.ed.ac.uk)
 *  (c) 2009 The University of Edinburgh
 *
 ****************************************************************************/

#ifndef STATS_FREE_ENERGY_H
#define STATS_FREE_ENERGY_H

#include "field.h"
#include "field_grad.h"
#include "map.h"

int stats_free_energy_density(field_t * q, map_t * map, int ncolloid);
int blue_phase_stats(field_t * q, field_grad_t * dq, map_t * map, int tstep);

#endif
