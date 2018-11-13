/*****************************************************************************
 *
 *  polar_active_rt.h
 *
 *  $Id: polar_active_rt.h 2479 2014-09-25 11:24:42Z stratford $
 *
 *  Edinburgh Soft Matter and Statistical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *  Kevin Stratford (kevin@epcc.ed.ac.uk)
 *  (c) 2011 The University of Edinburgh
 *
 *****************************************************************************/

#ifndef POLAR_ACTIVE_RT_H
#define POLAR_ACTIVE_RT_H

#include "field.h"

void polar_active_run_time(void);
int polar_active_rt_initial_conditions(field_t * p);
int polar_active_init_aster(field_t * p);

#endif
