/*****************************************************************************
 *
 *  blue_phase_rt.h
 *
 *  $Id: blue_phase_rt.h 2479 2014-09-25 11:24:42Z stratford $
 *
 *  Edinburgh Soft Matter and Statistical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *  Kevin Stratford (kevin@epcc.ed.ac.uk)
 *  (c) 2011 The University of Edinbrugh
 *
 *****************************************************************************/

#ifndef BLUE_PHASE_RT_H
#define BLUE_PHASE_RT_H

#include "field.h"

int blue_phase_rt_initial_conditions(field_t * q);
void blue_phase_run_time(void);

#endif
