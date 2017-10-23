/*****************************************************************************
 *
 *  phi_lb_coupler.h
 *
 *  Edinburgh Soft Matter and Statistical Physics Group and
 *  Edinburgh Parallel Computing Centre
 *
 *  (c) 2009-2014 The University of Edinburgh
 *  Kevin Stratford (kevin@epcc.ed.ac.uk)
 *
 *****************************************************************************/

#ifndef PHI_LB_COUPLER_H
#define PHI_LB_COUPLER_H

#include "field.h"
#include "model.h"

__targetHost__ int phi_lb_to_field(field_t * phi, lb_t * lb);
__targetHost__ int phi_lb_to_field_host(field_t * phi, lb_t * lb);
__targetHost__ int phi_lb_from_field(field_t * phi, lb_t * lb);

#endif
