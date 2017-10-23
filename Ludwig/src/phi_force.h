/*****************************************************************************
 *
 *  phi_force.h
 *
 *  $Id: phi_force.h 2636 2015-04-27 14:30:38Z agray3 $
 *
 *  Edinburgh Soft Matter and Statistical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *  Kevin Stratford (kevin@epcc.ed.ac.uk)
 *  (c) 2011 The University of Edinburgh
 *
 *****************************************************************************/

#ifndef PHI_FORCE_H
#define PHI_FORCE_H

#include "field.h"
#include "hydro.h"
#include "field_grad_s.h"

__targetHost__ int phi_force_calculation(field_t * phi, field_t* q_, field_grad_t* q_grad_, hydro_t * hydro);
__targetHost__ int phi_force_required(int * flag);
__targetHost__ int phi_force_required_set(const int flag);
__targetHost__ int phi_force_divergence_set(const int flag);

#endif
